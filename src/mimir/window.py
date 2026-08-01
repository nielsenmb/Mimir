"""Spectral-window calculation for sampled astronomical time series."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from operator import index
from typing import Any

import astropy.units as u
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import simpson

from mimir.inputs import _resolve_timeseries_input
from mimir.timeseries import TimeSeries


@dataclass(frozen=True)
class SpectralWindow:
    """A normalized spectral window centred on zero frequency.

    Parameters
    ----------
    frequency : numpy.ndarray
        Frequency offsets on a regular grid centred on zero.
    power : numpy.ndarray
        Squared modulus of the normalized Fourier transform of the sampling
        window.
    frequency_spacing : float
        Separation between adjacent frequency bins.
    effective_frequency_spacing : float
        Integral of the spectral window over the returned frequency range.
    nominal_frequency_spacing : float
        Reciprocal of the time-series duration.
    frequency_unit : str
        Unit shared by all frequency quantities.
    oversampling : int
        Number of samples per nominal Fourier spacing.
    half_width : float
        Maximum absolute returned frequency.
    """

    frequency: NDArray[np.float64]
    power: NDArray[np.float64]
    frequency_spacing: float
    effective_frequency_spacing: float
    nominal_frequency_spacing: float
    frequency_unit: str
    oversampling: int
    half_width: float

    @property
    def n_bins(self) -> int:
        """Number of frequency bins."""
        return int(self.frequency.size)


def spectral_window(
    time: ArrayLike | None = None,
    *,
    time_series: TimeSeries | None = None,
    target: str | None = None,
    mast_kwargs: Mapping[str, Any] | None = None,
    half_width: float | None = None,
    oversampling: int = 10,
    time_unit: str = "d",
    frequency_unit: str = "uHz",
) -> SpectralWindow:
    """Compute the normalized spectral window of a sampling pattern.

    Parameters
    ----------
    time : array-like, optional
        Sample times. Values are validated, sorted, and interpreted using
        ``time_unit``. Flux values are unnecessary for a sampling window.
    time_series : TimeSeries, optional
        Existing validated time series. Values of other types are rejected.
    target : str, optional
        Target name or identifier understood by Lightkurve. A target triggers
        MAST retrieval.
    mast_kwargs : mapping, optional
        Options passed to :func:`mimir.load_lightcurve` when ``target`` is
        selected. Put Lightkurve search constraints in the nested
        ``search_kwargs`` mapping.
    half_width : float, optional
        Maximum absolute frequency offset in ``frequency_unit``. The default
        is 100 times the nominal spacing, ``1 / T``.
    oversampling : int, default=10
        Number of samples per nominal Fourier spacing.
    time_unit : str, default="d"
        Astropy-compatible unit for array-like sample times.
    frequency_unit : str, default="uHz"
        Astropy-compatible unit for returned frequency quantities.

    Returns
    -------
    SpectralWindow
        Symmetric frequency offsets, normalized window power, and nominal and
        effective frequency spacings.

    Raises
    ------
    TypeError
        If exactly one of ``time``, ``time_series``, or ``target`` is not
        selected.

    Notes
    -----
    The window is the squared modulus of the discrete Fourier transform of
    unit weights at the observation times. Its zero-frequency value is one.
    The effective spacing is the numerical integral of the returned window.
    """
    series = _resolve_timeseries_input(
        time=time,
        flux=None,
        flux_err=None,
        time_series=time_series,
        target=target,
        mast_kwargs=mast_kwargs,
        time_unit=time_unit,
        flux_unit=None,
        allow_time_only=True,
    )
    oversampling_value = _validate_oversampling(oversampling)
    scale, frequency_unit_label = _frequency_conversion(
        series.time_unit,
        frequency_unit,
    )

    nominal_native = 1.0 / series.duration
    nominal_spacing = nominal_native * scale
    bin_spacing = nominal_spacing / oversampling_value
    if half_width is None:
        half_width_value = 100.0 * nominal_spacing
    else:
        half_width_value = _positive_finite("half_width", half_width)

    bins_each_side = int(np.floor(half_width_value / bin_spacing))
    if bins_each_side < 1:
        raise ValueError(
            "half_width must contain at least one frequency bin on each side "
            "of zero"
        )
    frequency = (
        np.arange(-bins_each_side, bins_each_side + 1, dtype=float) * bin_spacing
    )
    frequency_native = frequency / scale
    power = _sampling_window_power(series.time, frequency_native)
    effective_spacing = float(simpson(power, x=frequency))

    return SpectralWindow(
        frequency=frequency,
        power=power,
        frequency_spacing=float(bin_spacing),
        effective_frequency_spacing=effective_spacing,
        nominal_frequency_spacing=float(nominal_spacing),
        frequency_unit=frequency_unit_label,
        oversampling=oversampling_value,
        half_width=float(frequency[-1]),
    )


def effective_frequency_spacing(
    time: ArrayLike | None = None,
    *,
    time_series: TimeSeries | None = None,
    target: str | None = None,
    mast_kwargs: Mapping[str, Any] | None = None,
    half_width: float | None = None,
    oversampling: int = 10,
    time_unit: str = "d",
    frequency_unit: str = "uHz",
) -> float:
    """Estimate independent frequency-bin spacing from the spectral window.

    Parameters
    ----------
    time : array-like, optional
        Sample times.
    time_series : TimeSeries, optional
        Existing validated time series.
    target : str, optional
        Target name or identifier understood by Lightkurve.
    mast_kwargs : mapping, optional
        Options passed to :func:`mimir.load_lightcurve` when ``target`` is
        selected.
    half_width : float, optional
        Maximum absolute integration frequency in ``frequency_unit``.
    oversampling : int, default=10
        Number of samples per nominal Fourier spacing.
    time_unit : str, default="d"
        Astropy-compatible unit for array-like sample times.
    frequency_unit : str, default="uHz"
        Astropy-compatible unit for the returned spacing.

    Returns
    -------
    float
        Integral of the normalized spectral window in ``frequency_unit``.
    """
    return spectral_window(
        time=time,
        time_series=time_series,
        target=target,
        mast_kwargs=mast_kwargs,
        half_width=half_width,
        oversampling=oversampling,
        time_unit=time_unit,
        frequency_unit=frequency_unit,
    ).effective_frequency_spacing

def _sampling_window_power(
    time: NDArray[np.float64],
    frequency: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Evaluate the sampling-window power in bounded-memory chunks."""
    centred_time = time - np.mean(time)
    power = np.empty(frequency.size, dtype=float)
    chunk_size = max(1, min(256, int(2_000_000 / time.size)))
    for start in range(0, frequency.size, chunk_size):
        stop = min(start + chunk_size, frequency.size)
        phase = np.multiply.outer(frequency[start:stop], centred_time)
        transform = np.mean(np.exp(-2j * np.pi * phase), axis=1)
        power[start:stop] = np.square(np.abs(transform))
    return power


def _validate_oversampling(value: int) -> int:
    """Validate and return an integer oversampling factor."""
    if isinstance(value, bool):
        raise TypeError("oversampling must be a positive integer")
    try:
        result = index(value)
    except TypeError as error:
        raise TypeError("oversampling must be a positive integer") from error
    if result < 1:
        raise ValueError("oversampling must be a positive integer")
    return result


def _positive_finite(name: str, value: float) -> float:
    """Validate and return a positive finite scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a positive finite number") from error
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return result


def _frequency_conversion(time_unit: str, frequency_unit: str) -> tuple[float, str]:
    """Return the conversion from inverse time units to a frequency unit."""
    try:
        input_unit = u.Unit(time_unit) ** -1
    except (TypeError, ValueError) as error:
        raise ValueError(
            "time_unit must be an Astropy-compatible time unit"
        ) from error
    try:
        output_unit = u.Unit(frequency_unit)
        scale = input_unit.to(output_unit)
    except (TypeError, ValueError, u.UnitConversionError) as error:
        raise ValueError(
            "frequency_unit must be an Astropy-compatible frequency unit"
        ) from error
    return float(scale), output_unit.to_string()

