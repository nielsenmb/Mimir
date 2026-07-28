"""Characterisation tests for behaviour inherited from :mod:`pbjam.IO`.

These tests intentionally describe PBjam 2.0.4, including known defects. They
are a migration contract, not a specification of Mimir's final API.
"""

from types import SimpleNamespace

import numpy as np
import pytest

IO = pytest.importorskip("pbjam.IO")

pytestmark = pytest.mark.compatibility


def test_time_series_filters_invalid_and_user_masked_samples():
    """Invalid and explicitly masked samples should be removed together."""
    time = np.arange(6.0)
    flux = np.array([1.0, 2.0, np.nan, 4.0, 5.0, 6.0])
    flux_err = np.ones(6)
    bad_idx = np.array([False, False, False, True, False, False])

    result = IO.timeSeries(
        "target",
        lk_kwargs={},
        time=time,
        flux=flux,
        flux_err=flux_err,
        badIdx=bad_idx,
    )

    assert np.array_equal(result.time, [0.0, 1.0, 4.0, 5.0])
    assert np.array_equal(result.flux, [1.0, 2.0, 5.0, 6.0])
    assert np.array_equal(result.flux_err, np.ones(4))
    assert result.NT == 4
    assert result.dT == pytest.approx(5.0)
    assert result.dt == pytest.approx(1.0)
    assert result.dutyCycle == pytest.approx(0.8)


def test_regular_series_documents_legacy_duty_cycle_above_one():
    """PBjam's endpoint convention gives N/(N-1) for a complete series."""
    sample_count = 8
    result = IO.timeSeries(
        "target",
        lk_kwargs={},
        time=np.arange(sample_count, dtype=float),
        flux=np.ones(sample_count),
    )

    assert result.dutyCycle == pytest.approx(sample_count / (sample_count - 1))
    assert result.dutyCycle > 1.0


def test_unweighted_normalization_conforms_to_parseval():
    """Integrated power and power density should reproduce the variance."""
    cadence_days = 120.0 / 86400.0
    time = np.arange(512) * cadence_days
    frequency_hz = 800e-6
    flux = 20.0 * np.sin(2.0 * np.pi * frequency_hz * time * 86400.0)

    result = IO.psd("target", time=time, flux=flux)
    result(oversampling=2, method="fast")

    variance = np.mean((flux - np.mean(flux)) ** 2)
    frequency = np.asarray(result.freq)
    power = np.asarray(result.power)
    power_density = np.asarray(result.powerdensity)

    expected_spacing = np.full(frequency.size - 1, result.df * 1e6 / 2.0)
    assert np.diff(frequency) == pytest.approx(expected_spacing)
    assert frequency[-1] < result.Nyquist * 1e6
    assert np.sum(power) == pytest.approx(2.0 * variance, rel=2e-5)
    assert np.sum(power_density) * result.df * 1e6 == pytest.approx(
        variance,
        rel=2e-5,
    )
    assert frequency[np.argmax(power)] == pytest.approx(
        frequency_hz * 1e6,
        abs=result.df * 1e6,
    )


def test_amplitude_documents_legacy_formula():
    """Lock down the current dimensionally inconsistent amplitude formula."""
    cadence_days = 120.0 / 86400.0
    time = np.arange(512) * cadence_days
    flux = np.sin(2.0 * np.pi * 600e-6 * time * 86400.0)

    result = IO.psd("target", time=time, flux=flux)
    result(method="fast")

    power = np.asarray(result.power)
    expected = power * np.sqrt(power) / (2.0 * result.normfactor)
    assert np.asarray(result.amplitude) == pytest.approx(expected, rel=2e-5)


def test_weighted_normalization_uses_unweighted_arithmetic_mean():
    """Lock down PBjam's present weighted normalization convention."""
    result = IO.psd.__new__(IO.psd)
    result.ls = SimpleNamespace(
        t=np.arange(3.0),
        y=np.array([1.0, 2.0, 4.0]),
        dy=np.array([1.0, 2.0, 1.0]),
    )
    raw_power = np.array([2.0, 1.0])

    result._getNorm(raw_power)

    centered = result.ls.y - np.mean(result.ls.y)
    mean_square = np.sum((centered / result.ls.dy) ** 2)
    mean_square /= np.sum((1.0 / result.ls.dy) ** 2)
    assert result.normfactor == pytest.approx(mean_square / np.sum(raw_power))


def test_window_function_fills_internal_gap_and_optional_padding():
    """The time-domain window should represent gaps and requested padding."""
    result = IO.psd.__new__(IO.psd)
    result.TS = SimpleNamespace(time=np.array([1.0, 2.0, 4.0]), dt=1.0)

    time, window = result.getTSWindowFunction(tmin=0.0, tmax=6.0)

    assert np.array_equal(time, [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0])
    assert np.array_equal(window, [0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0])


def test_lightkurve_download_and_stitch_calls_are_separable(monkeypatch):
    """Archive tests should use mocks rather than make live MAST requests."""
    calls = {}

    class Values:
        """Minimal Lightkurve column stand-in."""

        def __init__(self, value):
            """Store a NumPy representation of the column."""
            self.value = np.asarray(value)

    class DummyLightCurve:
        """Minimal stitched light curve."""

        time = Values([0.0, 1.0, 2.0])
        flux = Values([1.0, 0.9, 1.1])
        flux_err = Values([0.1, 0.1, 0.2])

    class DummyCollection:
        """Minimal downloaded light-curve collection."""

        def stitch(self):
            """Record stitching and return the combined light curve."""
            calls["stitched"] = True
            return DummyLightCurve()

    class DummySearch:
        """Minimal Lightkurve search result."""

        def download_all(self, download_dir=None):
            """Record the cache directory and return downloaded curves."""
            calls["download_dir"] = download_dir
            return DummyCollection()

    def search_lightcurve(target, **kwargs):
        """Record search arguments and return a mock result."""
        calls["target"] = target
        calls["kwargs"] = kwargs
        return DummySearch()

    result = IO.timeSeries.__new__(IO.timeSeries)
    result.ID = "KIC 123"
    result.lk_kwargs = {"mission": "Kepler"}
    result.downloadDir = "/tmp/lightkurve"
    result.cleanLC = lambda light_curve, threshold: light_curve

    monkeypatch.setattr(IO.time, "sleep", lambda _: None)
    monkeypatch.setattr(IO.np.random, "uniform", lambda low, high: 0.0)
    monkeypatch.setattr(IO.lk, "search_lightcurve", search_lightcurve)

    time, flux, flux_err = result._getTS(outlierRejection=4)

    assert calls == {
        "target": "KIC 123",
        "kwargs": {"mission": "Kepler"},
        "download_dir": "/tmp/lightkurve",
        "stitched": True,
    }
    assert np.allclose(time, [0.0, 1.0, 2.0])
    assert np.allclose(flux, [1.0, 0.9, 1.1])
    assert np.allclose(flux_err, [0.1, 0.1, 0.2])
