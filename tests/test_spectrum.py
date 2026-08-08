"""Tests for power-spectrum calculation."""

import numpy as np
import pytest

from mimir import PowerSpectrum, TimeSeries, inputs, power_spectrum


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


def test_explicit_frequency_grid_is_returned_exactly_for_different_targets():
    """A shared explicit grid should be identical for different time series."""
    time_short, flux_short = _sine_series(sample_count=256)
    time_long, flux_long = _sine_series(sample_count=640)
    frequency = np.arange(100.0, 2000.0, 5.0)

    short = power_spectrum(time_short, flux_short, frequency=frequency)
    long = power_spectrum(time_long, flux_long, frequency=frequency)

    assert np.array_equal(short.frequency, frequency)
    assert np.array_equal(long.frequency, frequency)
    assert short.frequency_spacing == pytest.approx(5.0)
    assert long.frequency_spacing == pytest.approx(5.0)
    assert short.power.shape == frequency.shape
    assert long.power.shape == frequency.shape


def test_explicit_frequency_grid_recovers_injected_signal():
    """An explicit grid should evaluate the spectrum at its requested bins."""
    time, flux = _sine_series()
    frequency = np.arange(100.0, 2000.0, 2.0)

    result = power_spectrum(time, flux, frequency=frequency)

    assert result.frequency[np.argmax(result.power)] == pytest.approx(800.0)
    assert result.nyquist_factor == pytest.approx(
        frequency[-1] / result.nyquist_frequency
    )
    assert result.oversampling == pytest.approx(
        1.0 / (time[-1] * result.frequency_spacing / 1e6 * 86400.0)
    )


@pytest.mark.parametrize(
    ("frequency", "message"),
    [
        ([100.0], "at least two"),
        ([[100.0, 200.0]], "one-dimensional"),
        ([0.0, 100.0], "positive and finite"),
        ([100.0, np.inf], "positive and finite"),
        ([200.0, 100.0], "strictly increasing"),
        ([100.0, 200.0, 350.0], "regularly spaced"),
    ],
)
def test_invalid_explicit_frequency_grid_raises(frequency, message):
    """Explicit grids unsupported by nifty-ls should fail descriptively."""
    time, flux = _sine_series()

    with pytest.raises(ValueError, match=message):
        power_spectrum(time, flux, frequency=frequency)


@pytest.mark.parametrize("kwargs", [{"oversampling": 2}, {"nyquist_factor": 1.5}])
def test_explicit_frequency_grid_rejects_automatic_controls(kwargs):
    """Explicit and automatic frequency-grid controls should not be mixed."""
    time, flux = _sine_series()

    with pytest.raises(ValueError, match="cannot be changed"):
        power_spectrum(time, flux, frequency=[100.0, 200.0], **kwargs)


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

    result = power_spectrum(time_series=series)

    assert result.flux_unit == "ppm"
    with pytest.raises(TypeError, match="must be omitted"):
        power_spectrum(flux=flux, time_series=series)


def test_target_name_can_be_used_directly(monkeypatch):
    """A resolvable target should be accepted by the numerical entry point."""
    time, flux = _sine_series()
    expected = TimeSeries(time, flux, time_unit="d", flux_unit="ppm")
    calls = {}

    def load(target, kwargs):
        """Record the simulated MAST request."""
        calls.update(target=target, kwargs=kwargs)
        return expected

    monkeypatch.setattr(inputs, "_load_mast_target", load)

    result = power_spectrum(
        target="KIC 8006161",
        mast_kwargs={"search_kwargs": {"mission": "Kepler", "exptime": 60}},
    )

    assert result.flux_unit == "ppm"
    assert calls == {
        "target": "KIC 8006161",
        "kwargs": {"search_kwargs": {"mission": "Kepler", "exptime": 60}},
    }


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


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"time": [0.0, 1.0], "flux": [1.0, 2.0], "target": "KIC 1"},
        {
            "time_series": TimeSeries([0.0, 1.0], [1.0, 2.0]),
            "target": "KIC 1",
        },
    ],
)
def test_exactly_one_spectrum_input_form_is_required(kwargs):
    """The spectrum API should reject missing and ambiguous input forms."""
    with pytest.raises(TypeError, match="exactly one input form"):
        power_spectrum(**kwargs)


def test_time_series_and_target_types_are_strict():
    """Named input forms should reject values of the wrong semantic type."""
    with pytest.raises(TypeError, match="time_series must be"):
        power_spectrum(time_series=[0.0, 1.0])
    with pytest.raises(TypeError, match="target must be"):
        power_spectrum(target=123)


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


def test_grid_requires_a_physical_sub_nyquist_bin():
    """A super-Nyquist-only grid cannot define the Parseval scale."""
    with pytest.raises(ValueError, match="physical one-sided"):
        power_spectrum(
            [0.0, 1.0],
            [0.0, 1.0],
            nyquist_factor=4.0,
            frequency_unit="1/d",
        )
