"""Tests for shared time-series input resolution."""

import numpy as np
import pytest

from mimir import TimeSeries, as_timeseries, inputs


def test_existing_timeseries_is_reused():
    """An existing validated object should pass through unchanged."""
    expected = TimeSeries([0.0, 1.0], [1.0, 2.0])

    result = as_timeseries(time_series=expected)

    assert result is expected


def test_array_input_constructs_timeseries():
    """Separate arrays should be converted into a validated time series."""
    result = as_timeseries(
        time=[0.0, 1.0, 2.0],
        flux=[1.0, 0.0, 1.0],
        flux_err=[0.1, 0.2, 0.1],
        flux_unit="ppm",
    )

    assert np.array_equal(result.time, [0.0, 1.0, 2.0])
    assert np.array_equal(result.flux, [1.0, 0.0, 1.0])
    assert np.array_equal(result.flux_err, [0.1, 0.2, 0.1])
    assert result.flux_unit == "ppm"


def test_target_input_loads_from_mast(monkeypatch):
    """A target identifier should invoke the existing high-level MAST loader."""
    expected = TimeSeries([0.0, 1.0], [1.0, 2.0])
    calls = {}

    def load(target, kwargs):
        """Record the target and loader options."""
        calls.update(target=target, kwargs=kwargs)
        return expected

    monkeypatch.setattr(inputs, "_load_mast_target", load)

    result = as_timeseries(
        target="KIC 8006161",
        mast_kwargs={
            "search_kwargs": {"mission": "Kepler", "exptime": 60},
            "numax": 3500,
        },
    )

    assert result is expected
    assert calls == {
        "target": "KIC 8006161",
        "kwargs": {
            "search_kwargs": {"mission": "Kepler", "exptime": 60},
            "numax": 3500,
        },
    }


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {
            "time": [0.0, 1.0],
            "flux": [1.0, 2.0],
            "time_series": TimeSeries([0.0, 1.0], [1.0, 2.0]),
        },
        {"time": [0.0, 1.0], "flux": [1.0, 2.0], "target": "KIC 1"},
        {
            "time_series": TimeSeries([0.0, 1.0], [1.0, 2.0]),
            "target": "KIC 1",
        },
    ],
)
def test_exactly_one_input_form_is_required(kwargs):
    """Missing or ambiguous input forms should be rejected."""
    with pytest.raises(TypeError, match="exactly one input form"):
        as_timeseries(**kwargs)


def test_named_input_types_are_strict():
    """Object and target inputs should accept only their documented types."""
    with pytest.raises(TypeError, match="time_series must be"):
        as_timeseries(time_series=[0.0, 1.0])
    with pytest.raises(TypeError, match="target must be"):
        as_timeseries(target=123)
    with pytest.raises(ValueError, match="must not be empty"):
        as_timeseries(target="  ")


def test_form_specific_arguments_cannot_be_mixed():
    """Arguments belonging to another input form should be rejected."""
    series = TimeSeries([0.0, 1.0], [1.0, 2.0])
    with pytest.raises(TypeError, match="must be omitted"):
        as_timeseries(flux=[1.0, 2.0], time_series=series)
    with pytest.raises(TypeError, match="must be omitted"):
        as_timeseries(time_series=series, mast_kwargs={})
    with pytest.raises(TypeError, match="must be omitted"):
        as_timeseries(flux=[1.0, 2.0], target="KIC 1")
    with pytest.raises(TypeError, match="only be used with target"):
        as_timeseries(time=[0.0, 1.0], flux=[1.0, 2.0], mast_kwargs={})


def test_array_input_requires_flux():
    """Raw sample times require corresponding flux values."""
    with pytest.raises(TypeError, match="flux is required"):
        as_timeseries(time=[0.0, 1.0])
