"""Independent regression checks for normalization, masks, and time conversion."""

from types import SimpleNamespace
from unittest.mock import Mock

import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time, TimeDelta
from astropy.timeseries import LombScargle
from astropy.utils.masked import Masked

from mimir import (
    TimeSeries,
    lightcurve_to_timeseries,
    mast,
    power_spectrum,
    spectral_window,
)
from mimir import spectrum as spectrum_module


@pytest.mark.parametrize("weighted", [False, True])
def test_shared_peak_density_is_independent_of_requested_spacing(weighted):
    """Compare the review reproducer with an independent direct LS calculation.

    Parameters
    ----------
    weighted : bool
        Whether to use heteroscedastic uncertainties.
    """
    time = np.arange(2048, dtype=float)
    flux = 20 * np.sin(2 * np.pi * 0.12345 * time)
    error = np.linspace(0.5, 2.0, time.size) if weighted else None
    weights = np.ones(time.size) if error is None else error**-2
    residual = flux - np.average(flux, weights=weights)
    variance = np.average(residual**2, weights=weights)
    reference_spacing = 1 / np.ptp(time)
    reference_frequency = np.arange(1, 1024) * reference_spacing
    ls = LombScargle(time, flux, error, normalization="psd")
    reference_power = ls.power(reference_frequency, method="cython")
    expected = ls.power(0.12345, method="cython") * variance
    expected /= np.sum(reference_power) * reference_spacing
    for step in [0.0001, 0.0005, 0.001, 0.005, 0.01]:
        frequency = 0.12345 + np.arange(20) * step
        result = power_spectrum(
            time, flux, error, time_unit="s", frequency_unit="Hz",
            frequency=frequency, nthreads=1,
        )
        assert result.power_density[0] == pytest.approx(expected, rel=1e-7)
        np.testing.assert_array_equal(result.frequency, frequency)
        np.testing.assert_allclose(result.power, result.power_density * step)


def test_low_frequency_density_is_invariant_to_oversampling_and_band_limits():
    """Keep low-frequency density fixed when the output grid or band changes."""
    time = np.arange(1024, dtype=float)
    flux = np.sin(2 * np.pi * 1.25 * time / time[-1])
    base = power_spectrum(time, flux, time_unit="s", frequency_unit="Hz", nthreads=1)
    for factor in [0.5, 1.0, 1.5]:
        dense = power_spectrum(
            time, flux, time_unit="s", frequency_unit="Hz", oversampling=4,
            nyquist_factor=factor, nthreads=1,
        )
        shared = min(base.n_bins, dense.n_bins // 4)
        np.testing.assert_allclose(
            dense.power_density[3:shared * 4:4], base.power_density[:shared],
            rtol=1e-6, atol=1e-8,
        )


def test_dense_grid_forwards_threads_to_fixed_reference(monkeypatch):
    """Ensure the added reference evaluation receives the requested thread count.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Backend replacement fixture.
    """
    backend = Mock(side_effect=lambda *args, **kw: SimpleNamespace(
        power=np.ones(kw["Nf"]), backend="finufft",
    ))
    monkeypatch.setattr(spectrum_module.nifty_ls, "lombscargle", backend)
    power_spectrum(np.arange(64), np.sin(np.arange(64)), oversampling=4, nthreads=2)
    assert len(backend.call_args_list) == 2
    requested, reference = [call.kwargs for call in backend.call_args_list]
    assert requested["nthreads"] == reference["nthreads"] == 2
    assert reference["fmin"] == pytest.approx(1 / 63)
    assert requested["fmin"] == pytest.approx(reference["fmin"] / 4)


@pytest.mark.parametrize("array_type", [np.ma.array, Masked])
def test_masks_combine_with_bad_mask_before_strict_validation(array_type):
    """Drop masked invalid samples and retain correctly paired sorted measurements.

    Parameters
    ----------
    array_type : callable
        NumPy or Astropy masked-array constructor.
    """
    time = array_type([np.nan, 4., 3., 2., 1., 0.], mask=[1, 0, 0, 0, 0, 0])
    flux = array_type([0., np.nan, 2., 3., 4., 5.], mask=[0, 1, 0, 0, 0, 0])
    error = array_type([1., 1., -1., 1., 1., 1.], mask=[0, 0, 1, 0, 0, 0])
    result = TimeSeries(
        time, flux, error, bad_mask=[False, False, False, True, False, False],
        invalid="raise",
    )
    np.testing.assert_array_equal(result.time, [0., 1.])
    np.testing.assert_array_equal(result.flux, [5., 4.])
    assert result.input_size == 6
    assert result.n_removed == 4


def test_time_only_window_preserves_input_mask():
    """Exclude masked timestamps when the window API supplies synthetic flux."""
    time = np.ma.array(np.arange(32.), mask=False)
    time.mask[10:20] = True
    masked = spectral_window(time, time_unit="s", frequency_unit="Hz")
    clean = spectral_window(time.compressed(), time_unit="s", frequency_unit="Hz")
    np.testing.assert_allclose(masked.power, clean.power)


def test_conversion_excludes_masked_columns_before_ppm_normalization():
    """Exclude masked time, flux, and errors before estimating the flux median."""
    time = Time(2451545 + np.arange(6.), format="jd", scale="tdb")
    time[0] = np.ma.masked
    flux = Masked([100., 1000., 100., 10., 11., 9.], mask=[0, 1, 0, 0, 0, 0])
    error = np.ma.array(np.ones(6), mask=[0, 0, 1, 0, 0, 0])
    lc = SimpleNamespace(time=time, flux=flux, flux_err=error)
    result = lightcurve_to_timeseries(lc)
    np.testing.assert_allclose(result.time, [3., 4., 5.])
    np.testing.assert_allclose(result.flux, [0., 100000., -100000.])
    assert result.n_removed == 3
    np.testing.assert_array_equal(flux.unmasked, [100., 1000., 100., 10., 11., 9.])


@pytest.mark.parametrize("time_format", ["jd", "mjd", "unix", "isot"])
def test_time_formats_and_output_units_give_identical_spectra(time_format):
    """Compare the same absolute observations in different formats and units.

    Parameters
    ----------
    time_format : str
        Astropy display format to test.
    """
    seconds = 120 * np.arange(256, dtype=float)
    time = Time(2459000., seconds / 86400, format="jd", scale="tdb")
    time.format = time_format
    flux = 1 + 1e-4 * np.sin(2 * np.pi * 0.0008 * seconds)
    lc = SimpleNamespace(time=time, flux=flux, flux_err=np.full(time.size, 1e-5))
    days = lightcurve_to_timeseries(lc)
    sec = lightcurve_to_timeseries(lc, time_unit="s")
    assert days.cadence * 86400 == pytest.approx(120, abs=1e-6)
    assert sec.cadence == pytest.approx(120, abs=1e-6)
    assert days.time[0] == pytest.approx(7455)
    assert time.format == time_format
    assert mast._infer_exposure_time(lc) == pytest.approx(120, abs=1e-6)
    # Use an explicit common output grid to avoid endpoint rounding changing bin counts.
    frequency = np.arange(100., 2000., 10.)
    a = power_spectrum(time_series=days, frequency=frequency, nthreads=1)
    b = power_spectrum(time_series=sec, frequency=frequency, nthreads=1)
    np.testing.assert_allclose(a.power_density, b.power_density, rtol=1e-5, atol=1e-5)


@pytest.mark.parametrize(
    "time", [np.array([0., 1., 2.]), np.array([0., 86400., 172800.]) * u.s,
             TimeDelta([0., 1., 2.], format="jd")],
)
def test_relative_time_columns_are_converted(time):
    """Apply real unit conversion to relative coordinates and unitless day values.

    Parameters
    ----------
    time : array-like, Quantity, or TimeDelta
        Equivalent relative time coordinates.
    """
    result = lightcurve_to_timeseries(SimpleNamespace(time=time, flux=[1., 2., 3.]),
                                     time_unit="s")
    np.testing.assert_allclose(result.time, [0., 86400., 172800.])
    with pytest.raises(ValueError, match="time unit"):
        lightcurve_to_timeseries(SimpleNamespace(time=time, flux=[1., 2., 3.]),
                                time_unit="m")


@pytest.mark.mast
@pytest.mark.parametrize("time_format", ["bkjd", "btjd", "unix"])
def test_real_lightkurve_masks_and_mission_time_formats(time_format):
    """Exercise actual Lightkurve columns and its registered mission time formats.

    Parameters
    ----------
    time_format : str
        Registered mission or Unix time representation.
    """
    lk = pytest.importorskip("lightkurve")
    time = Time([2459000., 2459001., 2459002.], format="jd", scale="tdb")
    time.format = time_format
    lc = lk.LightCurve(
        time=time, flux=Masked([1., 1000., 1.], mask=[0, 1, 0]),
        flux_err=[0.1, 0.1, 0.1],
    )
    result = lightcurve_to_timeseries(lc)
    np.testing.assert_allclose(result.time, [7455., 7457.])
    np.testing.assert_array_equal(result.flux, [0., 0.])
    assert result.n_removed == 1
