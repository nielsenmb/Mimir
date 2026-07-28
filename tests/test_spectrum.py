"""Tests for power-spectrum calculation."""

import numpy as np
import pytest

from mimir import PowerSpectrum, TimeSeries, power_spectrum


def _sine_series(
    *,
    frequency_hz=800e-6,
    amplitude=20.0,
    sample_count=512,
    cadence_seconds=120.0,
):
    """Return a regularly sampled sinusoidal time series in days."""
    time = np.arange(sample_count) * cadence_seconds / 86400.0
    flux = amplitude * np.sin(2.0 * np.pi * frequency_hz * time * 86400.0)
    return time, flux


def test_sine_peak_and_frequency_grid():
    """A sinusoid should peak at the injected frequency on the requested grid."""
    time, flux = _sine_series()

    result = power_spectrum(time, flux, oversampling=2)

    expected_spacing = 1.0 / (time[-1] * 2.0) / 86400.0 * 1e6
    assert isinstance(result, PowerSpectrum)
    assert result.frequency_spacing == pytest.approx(expected_spacing)
    assert np.diff(result.frequency) == pytest.approx(expected_spacing)
    assert result.frequency[0] == pytest.approx(expected_spacing)
    assert result.maximum_frequency <= result.nyquist_frequency
    assert result.frequency[np.argmax(result.power)] == pytest.approx(
        800.0,
        abs=expected_spacing,
    )
    assert result.frequency_unit == "uHz"
    assert result.backend != "astropy"


def test_power_and_density_satisfy_parseval_relation():
    """Integrated power and power density should reproduce flux variance."""
    time, flux = _sine_series()

    result = power_spectrum(time, flux, oversampling=3)
    variance = np.var(flux)

    assert np.sum(result.power) == pytest.approx(variance, rel=1e-12)
    assert np.sum(result.power_density) * result.frequency_spacing == (
        pytest.approx(variance, rel=1e-12)
    )
    assert result.amplitude == pytest.approx(
        np.sqrt(2.0 * result.oversampling * result.power)
    )


def test_weighted_normalization_uses_weighted_mean_and_variance():
    """Uncertainties should define the Parseval target through their weights."""
    time, flux = _sine_series(sample_count=128)
    flux_err = np.linspace(0.5, 2.0, time.size)
    weights = 1.0 / flux_err**2
    mean = np.average(flux, weights=weights)
    weighted_variance = np.average((flux - mean) ** 2, weights=weights)

    result = power_spectrum(time, flux, flux_err)

    assert np.sum(result.power) == pytest.approx(weighted_variance, rel=1e-12)


def test_oversampling_and_nyquist_factor_control_grid():
    """Grid density and upper limit should follow their public parameters."""
    time, flux = _sine_series()

    base = power_spectrum(time, flux)
    dense = power_spectrum(time, flux, oversampling=4, nyquist_factor=1.5)

    assert dense.frequency_spacing == pytest.approx(base.frequency_spacing / 4.0)
    assert dense.nyquist_frequency == pytest.approx(base.nyquist_frequency)
    assert dense.maximum_frequency <= 1.5 * dense.nyquist_frequency
    assert dense.maximum_frequency > dense.nyquist_frequency


def test_frequency_units_are_converted_from_time_units():
    """Equivalent time units should produce identical physical frequencies."""
    time_days, flux = _sine_series()

    in_days = power_spectrum(time_days, flux, frequency_unit="Hz")
    in_seconds = power_spectrum(
        time_days * 86400.0,
        flux,
        time_unit="s",
        frequency_unit="Hz",
    )

    assert in_days.frequency == pytest.approx(in_seconds.frequency)
    assert in_days.frequency_spacing == pytest.approx(
        in_seconds.frequency_spacing
    )
    assert in_days.frequency_unit == "Hz"


def test_validated_time_series_can_be_reused():
    """The numerical layer should accept an existing TimeSeries directly."""
    time, flux = _sine_series()
    series = TimeSeries(time, flux, time_unit="d", flux_unit="ppm")

    result = power_spectrum(series)

    assert result.flux_unit == "ppm"
    with pytest.raises(TypeError, match="must be omitted"):
        power_spectrum(series, flux)


def test_constant_flux_returns_zero_spectrum():
    """A zero-variance series should have zero power without NaNs."""
    time = np.arange(32.0)

    result = power_spectrum(time, np.ones(time.size), frequency_unit="1 / d")

    assert np.all(result.power == 0.0)
    assert np.all(result.power_density == 0.0)
    assert np.all(result.amplitude == 0.0)


@pytest.mark.parametrize(
    ("kwargs", "exception", "message"),
    [
        ({"oversampling": 0}, ValueError, "positive integer"),
        ({"oversampling": 1.5}, TypeError, "positive integer"),
        ({"nyquist_factor": 0}, ValueError, "positive finite"),
        ({"nyquist_factor": np.inf}, ValueError, "positive finite"),
        ({"frequency_unit": "m"}, ValueError, "frequency unit"),
    ],
)
def test_invalid_spectrum_parameters_raise(kwargs, exception, message):
    """Invalid grid and unit options should produce descriptive exceptions."""
    time, flux = _sine_series()

    with pytest.raises(exception, match=message):
        power_spectrum(time, flux, **kwargs)


def test_array_input_requires_flux():
    """Raw time samples should not be accepted without flux values."""
    with pytest.raises(TypeError, match="flux is required"):
        power_spectrum(np.arange(8.0))



def _injected_bin(result):
    """Return the grid index nearest the injected frequency."""
    return int(np.argmin(np.abs(result.frequency - 800.0)))


def test_amplitude_is_stable_under_oversampling():
    """Zero-padding should refine the grid without diluting amplitude."""
    time, flux = _sine_series()

    base = power_spectrum(time, flux)
    dense = power_spectrum(time, flux, oversampling=4)

    assert dense.amplitude[_injected_bin(dense)] == pytest.approx(
        base.amplitude[_injected_bin(base)],
        rel=2e-3,
    )


def test_super_nyquist_bins_do_not_change_physical_band_normalization():
    """Aliased bins above Nyquist should not dilute the physical spectrum."""
    time, flux = _sine_series()

    base = power_spectrum(time, flux, oversampling=4)
    extended = power_spectrum(
        time,
        flux,
        oversampling=4,
        nyquist_factor=1.5,
    )
    base_index = _injected_bin(base)
    extended_index = _injected_bin(extended)

    assert extended.power_density[extended_index] == pytest.approx(
        base.power_density[base_index],
        rel=1e-9,
    )
    assert extended.amplitude[extended_index] == pytest.approx(
        base.amplitude[base_index],
        rel=1e-9,
    )
    assert np.sum(extended.power) > np.var(flux)


def test_under_resolved_frequency_grid_raises_descriptive_error():
    """A grid unsupported by nifty-ls should fail before backend dispatch."""
    with pytest.raises(ValueError, match="fewer than two bins"):
        power_spectrum(
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 0.0],
            frequency_unit="1/d",
        )
