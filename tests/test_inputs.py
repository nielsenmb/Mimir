"""Tests for shared time-series input coercion."""

import numpy as np
import pytest

from mimir import TimeSeries, as_timeseries, inputs


def test_existing_timeseries_is_reused():
    """An existing validated object should pass through unchanged."""
    expected = TimeSeries([0.0, 1.0], [1.0, 2.0])

    result = as_timeseries(expected)

    assert result is expected


def test_array_input_constructs_timeseries():
    """Separate arrays should be converted into a validated time series."""
    result = as_timeseries(
        [0.0, 1.0, 2.0],
        [1.0, 0.0, 1.0],
        [0.1, 0.2, 0.1],
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
        "KIC 8006161",
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
    ("source", "flux", "mast_kwargs", "message"),
    [
        (TimeSeries([0.0, 1.0], [1.0, 2.0]), [1.0, 2.0], None, "must be omitted"),
        (TimeSeries([0.0, 1.0], [1.0, 2.0]), None, {}, "must be omitted"),
        ("KIC 1", [1.0, 2.0], None, "must be omitted"),
        ([0.0, 1.0], None, None, "flux is required"),
        ([0.0, 1.0], [1.0, 2.0], {}, "target identifier"),
    ],
)
def test_incompatible_input_arguments_raise(source, flux, mast_kwargs, message):
    """Arguments from different input forms should not be mixed silently."""
    with pytest.raises(TypeError, match=message):
        as_timeseries(source, flux, mast_kwargs=mast_kwargs)
