"""Power-spectrum and power-density-spectrum calculation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from operator import index
from typing import Any

import astropy.units as u
import nifty_ls
import numpy as np
from numpy.typing import ArrayLike, NDArray

from mimir.inputs import as_timeseries
from mimir.timeseries import TimeSeries


@dataclass(frozen=True)
class PowerSpectrum:
    """A one-sided, Parseval-normalized power spectrum.

    Parameters
    ----------
    frequency : numpy.ndarray
        Positive frequencies on a regular grid.
    power : numpy.ndarray
        Power in each frequency bin, normalized against the physical
        one-sided band through Nyquist.
    power_density : numpy.ndarray
        Power per unit frequency.
    amplitude : numpy.ndarray
        Sinusoidal semi-amplitude, defined as
        ``sqrt(2 * oversampling * power)``.
    frequency_spacing : float
        Separation between adjacent frequency bins.
    nyquist_frequency : float
        Nyquist frequency estimated from the median cadence.
    frequency_unit : str
        Unit label shared by the frequency quantities.
    flux_unit : str, optional
        Unit label for amplitude. Power and power density carry the
        corresponding squared flux units.
    oversampling : float
        Effective number of frequency samples per nominal Fourier spacing.
        This is an integer for automatically generated grids and may be a
        non-integer for a user-supplied grid.
    nyquist_factor : float
        Requested maximum frequency as a multiple of the Nyquist frequency.
    backend : str
        Numerical backend selected by nifty-ls.

    Notes
    -----
    Oversampled periodogram bins are correlated. The oversampling factor in the
    amplitude conversion compensates for the narrower bins, so refining the
    frequency grid does not dilute a coherent sinusoid's semi-amplitude.
    """

    frequency: NDArray[np.float64]
    power: NDArray[np.float64]
    power_density: NDArray[np.float64]
    amplitude: NDArray[np.float64]
    frequency_spacing: float
    nyquist_frequency: float
    frequency_unit: str
    flux_unit: str | None
    oversampling: float
    nyquist_factor: float
    backend: str

    @property
    def n_bins(self) -> int:
        """Number of frequency bins."""
        return int(self.frequency.size)

    @property
    def maximum_frequency(self) -> float:
        """Largest evaluated frequency."""
        return float(self.frequency[-1])


def power_spectrum(
    time: ArrayLike | None = None,
    flux: ArrayLike | None = None,
    flux_err: ArrayLike | None = None,
    *,
    time_series: TimeSeries | None = None,
    target: str | None = None,
    mast_kwargs: Mapping[str, Any] | None = None,
    frequency: ArrayLike | None = None,
    oversampling: int = 1,
    nyquist_factor: float = 1.0,
    time_unit: str = "d",
    flux_unit: str | None = None,
    frequency_unit: str = "uHz",
    backend: str = "auto",
) -> PowerSpectrum:
    """Compute a one-sided power spectrum with nifty-ls.

    Parameters
    ----------
    time : array-like, optional
        Sample times. ``flux`` is required when this input form is selected.
    flux : array-like, optional
        Flux measurements corresponding to ``time``.
    flux_err : array-like, optional
        Positive one-sigma flux uncertainties corresponding to ``time``.
    time_series : TimeSeries, optional
        Existing validated time series. Values of other types are rejected.
    target : str, optional
        Target name or identifier understood by Lightkurve. A target triggers
        MAST retrieval.
    mast_kwargs : mapping, optional
        Options passed to :func:`mimir.load_lightcurve` when ``target`` is
        selected. Put Lightkurve search constraints in the nested
        ``search_kwargs`` mapping.
    frequency : array-like, optional
        Explicit positive, strictly increasing, regularly spaced frequency
        grid in ``frequency_unit``. The returned frequencies match these
        values exactly. When supplied, ``oversampling`` and
        ``nyquist_factor`` must retain their default values.
    oversampling : int, default=1
        Number of frequency samples per nominal Fourier spacing, ``1 / T``.
    nyquist_factor : float, default=1.0
        Maximum evaluated frequency as a multiple of the Nyquist frequency
        estimated from the median cadence.
    time_unit : str, default="d"
        Astropy-compatible unit of array-like sample times.
    flux_unit : str, optional
        Unit label for array-like flux values.
    frequency_unit : str, default="uHz"
        Astropy-compatible unit for returned frequencies.
    backend : str, default="auto"
        Backend passed to :func:`nifty_ls.lombscargle`.

    Returns
    -------
    PowerSpectrum
        Regular frequencies, power, power density, amplitude, and grid
        metadata.

    Raises
    ------
    TypeError
        If exactly one of array input, ``time_series``, or ``target`` is not
        selected.
    ValueError
        If ``frequency`` is invalid or is combined with automatic-grid
        controls.

    Notes
    -----
    Power is rescaled using the physical one-sided band through the Nyquist
    frequency. At ``nyquist_factor=1``, its sum satisfies Parseval's relation.
    Truncated spectra contain only the variance represented in the returned
    band, while super-Nyquist aliases do not dilute the physical spectrum.
    With uncertainties, inverse-variance weights define the target variance.
    """
    series = as_timeseries(
        time=time,
        flux=flux,
        flux_err=flux_err,
        time_series=time_series,
        target=target,
        mast_kwargs=mast_kwargs,
        time_unit=time_unit,
        flux_unit=flux_unit,
    )
    frequency_scale, frequency_unit_label = _frequency_conversion(
        series.time_unit,
        frequency_unit,
    )
    nyquist_native = 0.5 / series.cadence
    if frequency is None:
        (
            output_frequency,
            spacing_native,
            oversampling_value,
            nyquist_factor_value,
            raw_power,
            normalization_power,
            selected_backend,
        ) = _automatic_grid_power(
            series,
            oversampling=oversampling,
            nyquist_factor=nyquist_factor,
            frequency_scale=frequency_scale,
            nyquist_native=nyquist_native,
            backend=backend,
        )
    else:
        if oversampling != 1 or nyquist_factor != 1.0:
            raise ValueError(
                "oversampling and nyquist_factor cannot be changed when "
                "frequency is supplied"
            )
        output_frequency, spacing_output = _validate_frequency_grid(frequency)
        spacing_native = spacing_output / frequency_scale
        requested_native = output_frequency / frequency_scale
        oversampling_value = 1.0 / (series.duration * spacing_native)
        nyquist_factor_value = float(requested_native[-1] / nyquist_native)
        raw_power, selected_backend = _nifty_power(
            series,
            minimum_native=float(requested_native[0]),
            spacing_native=spacing_native,
            n_bins=output_frequency.size,
            backend=backend,
        )
        normalization_bins = _normalization_bin_count(
            nyquist_native,
            spacing_native,
        )
        normalization_power, _ = _nifty_power(
            series,
            minimum_native=spacing_native,
            spacing_native=spacing_native,
            n_bins=max(2, normalization_bins),
            backend=backend,
        )
        normalization_power = normalization_power[:normalization_bins]

    target_variance = _flux_variance(series)
    power = _parseval_power(
        raw_power,
        target_variance,
        normalization_power=normalization_power,
    )

    frequency_spacing = spacing_native * frequency_scale
    power_density = power / frequency_spacing
    amplitude = np.sqrt(2.0 * oversampling_value * power)

    return PowerSpectrum(
        frequency=output_frequency,
        power=power,
        power_density=power_density,
        amplitude=amplitude,
        frequency_spacing=float(frequency_spacing),
        nyquist_frequency=float(nyquist_native * frequency_scale),
        frequency_unit=frequency_unit_label,
        flux_unit=series.flux_unit,
        oversampling=oversampling_value,
        nyquist_factor=nyquist_factor_value,
        backend=selected_backend,
    )


def _automatic_grid_power(
    series: TimeSeries,
    *,
    oversampling: int,
    nyquist_factor: float,
    frequency_scale: float,
    nyquist_native: float,
    backend: str,
) -> tuple[
    NDArray[np.float64],
    float,
    int,
    float,
    NDArray[np.float64],
    NDArray[np.float64],
    str,
]:
    """Evaluate a spectrum on Mimir's automatically generated grid."""
    oversampling_value = _validate_oversampling(oversampling)
    nyquist_factor_value = _positive_finite("nyquist_factor", nyquist_factor)
    spacing_native = 1.0 / (series.duration * oversampling_value)
    maximum_native = nyquist_factor_value * nyquist_native
    n_bins = int(np.floor(maximum_native / spacing_native))
    if n_bins < 2:
        raise ValueError(
            "the requested frequency range contains fewer than two bins; "
            "increase oversampling or nyquist_factor, or provide a longer "
            "time series"
        )

    normalization_bins = _normalization_bin_count(nyquist_native, spacing_native)
    evaluated_bins = max(n_bins, normalization_bins)
    raw_power, selected_backend = _nifty_power(
        series,
        minimum_native=spacing_native,
        spacing_native=spacing_native,
        n_bins=evaluated_bins,
        backend=backend,
    )
    frequency_spacing = spacing_native * frequency_scale
    frequency = np.arange(1, n_bins + 1, dtype=float) * frequency_spacing
    return (
        frequency,
        spacing_native,
        oversampling_value,
        nyquist_factor_value,
        raw_power[:n_bins],
        raw_power[:normalization_bins],
        selected_backend,
    )


def _validate_frequency_grid(frequency: ArrayLike) -> tuple[NDArray[np.float64], float]:
    """Validate and copy a user-supplied regular frequency grid."""
    try:
        values = np.asarray(frequency, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("frequency must be a one-dimensional numeric array") from error
    if values.ndim != 1:
        raise ValueError("frequency must be one-dimensional")
    if values.size < 2:
        raise ValueError("frequency must contain at least two values")
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("frequency values must be positive and finite")
    differences = np.diff(values)
    if np.any(differences <= 0.0):
        raise ValueError("frequency values must be strictly increasing")
    spacing = float(differences[0])
    if not np.allclose(differences, spacing, rtol=1e-10, atol=0.0):
        raise ValueError("frequency must be regularly spaced for nifty-ls")
    return values.copy(), spacing


def _normalization_bin_count(nyquist_native: float, spacing_native: float) -> int:
    """Return the number of regular bins in the physical one-sided band."""
    normalization_bins = int(np.floor(nyquist_native / spacing_native))
    if normalization_bins < 1:
        raise ValueError(
            "the physical one-sided frequency band contains no bins; use a "
            "finer frequency grid or provide a longer time series"
        )
    return normalization_bins


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
    """Validate a positive finite scalar."""
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
            "time_unit must be an Astropy-compatible time unit for spectrum "
            "calculation"
        ) from error
    try:
        output_unit = u.Unit(frequency_unit)
        scale = input_unit.to(output_unit)
    except (TypeError, ValueError, u.UnitConversionError) as error:
        raise ValueError(
            "frequency_unit must be an Astropy-compatible frequency unit"
        ) from error
    return float(scale), output_unit.to_string()


def _nifty_power(
    series: TimeSeries,
    *,
    minimum_native: float,
    spacing_native: float,
    n_bins: int,
    backend: str,
) -> tuple[NDArray[np.float64], str]:
    """Evaluate nifty-ls on the requested regular frequency grid."""
    result = nifty_ls.lombscargle(
        series.time,
        series.flux,
        dy=series.flux_err,
        fmin=minimum_native,
        fmax=minimum_native + (n_bins - 1) * spacing_native,
        Nf=n_bins,
        center_data=True,
        fit_mean=True,
        normalization="psd",
        assume_sorted_t=True,
        backend=backend,
    )
    raw_power = np.asarray(result.power, dtype=float)
    if raw_power.shape != (n_bins,) or not np.all(np.isfinite(raw_power)):
        raise RuntimeError("nifty-ls returned an invalid power spectrum")
    return np.clip(raw_power, 0.0, None), str(result.backend)


def _flux_variance(series: TimeSeries) -> float:
    """Return the unweighted or inverse-variance-weighted flux variance."""
    if series.flux_err is None:
        return float(np.mean((series.flux - np.mean(series.flux)) ** 2))
    weights = 1.0 / np.square(series.flux_err)
    mean = np.average(series.flux, weights=weights)
    return float(np.average(np.square(series.flux - mean), weights=weights))


def _parseval_power(
    raw_power: NDArray[np.float64],
    target_variance: float,
    *,
    normalization_power: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Scale a periodogram using its physical one-sided frequency band.

    Parameters
    ----------
    raw_power : numpy.ndarray
        Non-negative periodogram values, possibly extending above Nyquist.
    target_variance : float
        Flux variance that the bins through Nyquist must reproduce.
    normalization_power : numpy.ndarray
        Periodogram values on a regular grid covering the physical one-sided
        band through Nyquist.

    Returns
    -------
    numpy.ndarray
        Scaled periodogram values.
    """
    if target_variance == 0.0:
        return np.zeros_like(raw_power)
    total = float(np.sum(normalization_power))
    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError("nifty-ls returned no positive finite power")
    return raw_power * (target_variance / total)
