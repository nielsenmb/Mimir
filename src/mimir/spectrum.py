"""Power-spectrum and power-density-spectrum calculation."""

from __future__ import annotations

from dataclasses import dataclass
from operator import index

import astropy.units as u
import nifty_ls
import numpy as np
from numpy.typing import ArrayLike, NDArray

from mimir.timeseries import TimeSeries


@dataclass(frozen=True)
class PowerSpectrum:
    """A one-sided, Parseval-normalized power spectrum.

    Parameters
    ----------
    frequency : numpy.ndarray
        Positive frequencies on a regular grid.
    power : numpy.ndarray
        Power in each frequency bin, normalized so that its sum equals the
        variance of the input flux.
    power_density : numpy.ndarray
        Power per unit frequency.
    amplitude : numpy.ndarray
        Sinusoidal semi-amplitude per bin, defined as
        ``sqrt(2 * power)``.
    frequency_spacing : float
        Separation between adjacent frequency bins.
    nyquist_frequency : float
        Nyquist frequency estimated from the median cadence.
    frequency_unit : str
        Unit label shared by the frequency quantities.
    flux_unit : str, optional
        Unit label for amplitude. Power and power density carry the
        corresponding squared flux units.
    oversampling : int
        Number of frequency samples per nominal Fourier spacing.
    nyquist_factor : float
        Requested maximum frequency as a multiple of the Nyquist frequency.
    backend : str
        Numerical backend selected by nifty-ls.

    Notes
    -----
    Oversampled periodogram bins are correlated. Consequently, ``amplitude``
    describes individual bins and should not be interpreted as a collection of
    independent Fourier amplitudes when ``oversampling`` is greater than one.
    """

    frequency: NDArray[np.float64]
    power: NDArray[np.float64]
    power_density: NDArray[np.float64]
    amplitude: NDArray[np.float64]
    frequency_spacing: float
    nyquist_frequency: float
    frequency_unit: str
    flux_unit: str | None
    oversampling: int
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
    time: TimeSeries | ArrayLike,
    flux: ArrayLike | None = None,
    flux_err: ArrayLike | None = None,
    *,
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
    time : TimeSeries or array-like
        A validated time series or the sample times. When an array is supplied,
        ``flux`` is required and the inputs are validated by
        :class:`mimir.TimeSeries`.
    flux : array-like, optional
        Flux measurements. Must be omitted when ``time`` is a ``TimeSeries``.
    flux_err : array-like, optional
        Positive one-sigma flux uncertainties. Must be omitted when ``time`` is
        a ``TimeSeries``.
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

    Notes
    -----
    Power is rescaled to satisfy Parseval's relation over the evaluated
    one-sided grid. With uncertainties, the target variance uses inverse
    variance weights and a weighted mean.
    """
    series = _as_time_series(time, flux, flux_err, time_unit, flux_unit)
    oversampling_value = _validate_oversampling(oversampling)
    nyquist_factor_value = _positive_finite("nyquist_factor", nyquist_factor)
    frequency_scale, frequency_unit_label = _frequency_conversion(
        series.time_unit,
        frequency_unit,
    )

    spacing_native = 1.0 / (series.duration * oversampling_value)
    nyquist_native = 0.5 / series.cadence
    maximum_native = nyquist_factor_value * nyquist_native
    n_bins = int(np.floor(maximum_native / spacing_native))
    if n_bins < 1:
        raise ValueError(
            "the requested frequency range contains no bins; increase "
            "nyquist_factor or provide a longer time series"
        )

    raw_power, selected_backend = _nifty_power(
        series,
        spacing_native=spacing_native,
        n_bins=n_bins,
        backend=backend,
    )
    target_variance = _flux_variance(series)
    power = _parseval_power(raw_power, target_variance)

    frequency_spacing = spacing_native * frequency_scale
    frequency = np.arange(1, n_bins + 1, dtype=float) * frequency_spacing
    power_density = power / frequency_spacing
    amplitude = np.sqrt(2.0 * power)

    return PowerSpectrum(
        frequency=frequency,
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


def _as_time_series(
    time: TimeSeries | ArrayLike,
    flux: ArrayLike | None,
    flux_err: ArrayLike | None,
    time_unit: str,
    flux_unit: str | None,
) -> TimeSeries:
    """Return a validated time series from either supported input form."""
    if isinstance(time, TimeSeries):
        if flux is not None or flux_err is not None:
            raise TypeError(
                "flux and flux_err must be omitted when time is a TimeSeries"
            )
        return time
    if flux is None:
        raise TypeError("flux is required when time is an array")
    return TimeSeries(
        time=time,
        flux=flux,
        flux_err=flux_err,
        time_unit=time_unit,
        flux_unit=flux_unit,
    )


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
    spacing_native: float,
    n_bins: int,
    backend: str,
) -> tuple[NDArray[np.float64], str]:
    """Evaluate nifty-ls on the requested regular frequency grid."""
    result = nifty_ls.lombscargle(
        series.time,
        series.flux,
        dy=series.flux_err,
        fmin=spacing_native,
        fmax=n_bins * spacing_native,
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
) -> NDArray[np.float64]:
    """Scale raw periodogram values so that their sum is the flux variance."""
    if target_variance == 0.0:
        return np.zeros_like(raw_power)
    total = float(np.sum(raw_power))
    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError("nifty-ls returned no positive finite power")
    return raw_power * (target_variance / total)
