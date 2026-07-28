"""Tests for time-series validation and metadata."""

import numpy as np
import pytest

from mimir import TimeSeries


def test_regular_series_metadata():
    """A complete regular series should have a unit duty cycle."""
    result = TimeSeries(
        time=np.arange(8.0),
        flux=np.linspace(0.9, 1.1, 8),
        flux_err=np.full(8, 0.01),
    )

    assert result.n_samples == 8
    assert result.input_size == 8
    assert result.n_removed == 0
    assert result.cadence == pytest.approx(1.0)
    assert result.duration == pytest.approx(7.0)
    assert result.duty_cycle == pytest.approx(1.0)
    assert result.has_uncertainties
    assert result.time_unit == "d"
    assert result.flux_unit is None


def test_invalid_and_masked_samples_are_removed_together():
    """Non-finite and explicitly bad samples should be removed."""
    result = TimeSeries(
        time=np.arange(6.0),
        flux=np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0]),
        flux_err=np.ones(6),
        bad_mask=np.array([False, False, False, True, False, False]),
    )

    assert np.array_equal(result.time, [0.0, 1.0, 4.0, 5.0])
    assert np.array_equal(result.flux, [1.0, 2.0, 5.0, 6.0])
    assert np.array_equal(result.flux_err, np.ones(4))
    assert result.input_size == 6
    assert result.n_removed == 2
    assert result.duty_cycle == pytest.approx(4.0 / 6.0)


def test_samples_are_sorted_with_all_columns_aligned():
    """Default sorting should preserve correspondence between columns."""
    result = TimeSeries(
        time=[2.0, 0.0, 1.0],
        flux=[20.0, 0.0, 10.0],
        flux_err=[0.2, 0.1, 0.15],
    )

    assert np.array_equal(result.time, [0.0, 1.0, 2.0])
    assert np.array_equal(result.flux, [0.0, 10.0, 20.0])
    assert np.array_equal(result.flux_err, [0.1, 0.15, 0.2])


def test_unordered_samples_can_be_rejected():
    """Callers should be able to require already ordered input."""
    with pytest.raises(ValueError, match="strictly increasing"):
        TimeSeries(time=[1.0, 0.0], flux=[1.0, 2.0], sort=False)


def test_unmasked_non_finite_samples_can_be_rejected():
    """Strict invalid-value handling should report non-finite samples."""
    with pytest.raises(ValueError, match="1 non-finite"):
        TimeSeries(time=[0.0, 1.0], flux=[1.0, np.nan], invalid="raise")


def test_masked_non_finite_samples_do_not_fail_strict_validation():
    """A caller-provided mask should take precedence over strict validation."""
    result = TimeSeries(
        time=[0.0, 1.0, 2.0],
        flux=[1.0, np.nan, 3.0],
        bad_mask=[False, True, False],
        invalid="raise",
    )

    assert np.array_equal(result.time, [0.0, 2.0])


@pytest.mark.parametrize(
    ("time", "flux", "message"),
    [
        ([[0.0, 1.0]], [1.0, 2.0], "time must be one-dimensional"),
        ([0.0, 1.0], [[1.0, 2.0]], "flux must be one-dimensional"),
        ([0.0, 1.0], [1.0], "same length"),
        ([0.0, 0.0], [1.0, 2.0], "duplicate"),
        ([np.nan, 1.0], [1.0, 2.0], "at least two valid"),
    ],
)
def test_invalid_core_inputs_raise(time, flux, message):
    """Malformed core inputs should produce descriptive exceptions."""
    with pytest.raises(ValueError, match=message):
        TimeSeries(time=time, flux=flux)


def test_flux_uncertainties_must_be_valid():
    """Uncertainties should match the data and be strictly positive."""
    with pytest.raises(ValueError, match="same length"):
        TimeSeries(time=[0.0, 1.0], flux=[1.0, 2.0], flux_err=[0.1])

    with pytest.raises(ValueError, match="strictly positive"):
        TimeSeries(time=[0.0, 1.0], flux=[1.0, 2.0], flux_err=[0.1, 0.0])


def test_bad_mask_must_be_boolean_and_match_data():
    """Bad masks should not silently coerce integers or change shape."""
    with pytest.raises(TypeError, match="boolean"):
        TimeSeries(time=[0.0, 1.0], flux=[1.0, 2.0], bad_mask=[0, 1])

    with pytest.raises(ValueError, match="same length"):
        TimeSeries(time=[0.0, 1.0], flux=[1.0, 2.0], bad_mask=[False])


def test_input_arrays_are_copied():
    """Later mutation of caller-owned arrays should not alter the result."""
    time = np.array([0.0, 1.0])
    flux = np.array([1.0, 2.0])
    result = TimeSeries(time=time, flux=flux)

    time[:] = 10.0
    flux[:] = 20.0

    assert np.array_equal(result.time, [0.0, 1.0])
    assert np.array_equal(result.flux, [1.0, 2.0])


def test_unit_labels_are_explicit_metadata():
    """Unit labels should describe arrays and derived metadata."""
    result = TimeSeries(
        time=[0.0, 120.0],
        flux=[1.0, 1.1],
        time_unit="s",
        flux_unit="relative flux",
    )

    assert result.time_unit == "s"
    assert result.flux_unit == "relative flux"
    assert result.cadence == pytest.approx(120.0)

    with pytest.raises(ValueError, match="time_unit"):
        TimeSeries(time=[0.0, 1.0], flux=[1.0, 2.0], time_unit="")
