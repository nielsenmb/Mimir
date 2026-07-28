"""Tests for optional Lightkurve and MAST integration."""

from types import SimpleNamespace

import numpy as np
import pytest

from mimir import mast


class Column:
    """Minimal Lightkurve column stand-in."""

    def __init__(self, values):
        """Store column values."""
        self.value = np.asarray(values)


class DummyLightCurve:
    """Record the reduction operations applied to a light curve."""

    def __init__(self):
        """Create a small regularly sampled light curve."""
        self.time = Column([0.0, 1.0, 2.0])
        self.flux = Column([1.0, 0.9, 1.1])
        self.flux_err = Column([0.1, 0.1, 0.2])
        self.meta = {"EXPTIME": 120.0}
        self.calls = []

    def normalize(self):
        """Record normalization."""
        self.calls.append(("normalize", {}))
        return self

    def remove_nans(self):
        """Record non-finite removal."""
        self.calls.append(("remove_nans", {}))
        return self

    def remove_outliers(self, **kwargs):
        """Record outlier removal."""
        self.calls.append(("remove_outliers", kwargs))
        return self

    def flatten(self, **kwargs):
        """Record flattening."""
        self.calls.append(("flatten", kwargs))
        return self


def test_search_lightcurves_forwards_arguments(monkeypatch):
    """Search arguments should be forwarded without modification."""
    calls = {}
    expected = [object()]

    def search_lightcurve(target, **kwargs):
        """Record the simulated Lightkurve search."""
        calls.update(target=target, kwargs=kwargs)
        return expected

    monkeypatch.setattr(
        mast,
        "_import_lightkurve",
        lambda: SimpleNamespace(search_lightcurve=search_lightcurve),
    )

    result = mast.search_lightcurves("KIC 123", mission="Kepler", exptime=60)

    assert result is expected
    assert calls == {
        "target": "KIC 123",
        "kwargs": {"mission": "Kepler", "exptime": 60},
    }


@pytest.mark.parametrize("target", ["", "  ", None])
def test_search_lightcurves_rejects_invalid_target(target):
    """Target identifiers should be non-empty strings."""
    with pytest.raises(ValueError, match="target"):
        mast.search_lightcurves(target)


def test_search_lightcurves_rejects_empty_result(monkeypatch):
    """An empty archive search should raise a clear lookup error."""
    monkeypatch.setattr(
        mast,
        "_import_lightkurve",
        lambda: SimpleNamespace(search_lightcurve=lambda target: []),
    )

    with pytest.raises(LookupError, match="no light curves"):
        mast.search_lightcurves("missing")


def test_import_lightkurve_explains_optional_extra(monkeypatch):
    """A missing optional dependency should name the installation extra."""

    def missing(name):
        """Simulate an unavailable Lightkurve package."""
        raise ImportError(name)

    monkeypatch.setattr(mast, "import_module", missing)

    with pytest.raises(ImportError, match=r"mimir-astro\[mast\]"):
        mast._import_lightkurve()


def test_download_lightcurves_forwards_cache_directory():
    """Downloads should preserve the caller's cache directory."""
    calls = {}
    collection = [object()]

    class Search:
        """Minimal non-empty search result."""

        def __len__(self):
            """Return one simulated product."""
            return 1

        def download_all(self, download_dir=None):
            """Record and return a simulated download."""
            calls["download_dir"] = download_dir
            return collection

    result = mast.download_lightcurves(Search(), download_dir="/cache")

    assert result is collection
    assert calls == {"download_dir": "/cache"}


@pytest.mark.parametrize("result", [None, []])
def test_download_lightcurves_rejects_empty_download(result):
    """Missing downloaded products should raise a lookup error."""

    class Search:
        """Minimal non-empty search result."""

        def __len__(self):
            """Return one simulated product."""
            return 1

        def download_all(self, download_dir=None):
            """Return the requested empty result."""
            return result

    with pytest.raises(LookupError, match="did not download"):
        mast.download_lightcurves(Search())


def test_reduce_lightcurve_stitches_and_applies_operations():
    """Collections should be stitched before the PBjam reduction sequence."""
    light_curve = DummyLightCurve()

    class Collection:
        """Minimal downloaded light-curve collection."""

        def stitch(self):
            """Return the simulated stitched light curve."""
            return light_curve

    result = mast.reduce_lightcurve(
        Collection(),
        outlier_sigma=4,
        flatten_window_length=20,
    )

    assert result is light_curve
    assert light_curve.calls == [
        ("normalize", {}),
        ("remove_nans", {}),
        ("remove_outliers", {"sigma": 4.0}),
        ("flatten", {"window_length": 21}),
    ]


def test_reduce_lightcurve_can_skip_optional_operations():
    """Normalization, outlier removal, and flattening should be configurable."""
    light_curve = DummyLightCurve()

    mast.reduce_lightcurve(
        light_curve,
        normalize=False,
        outlier_sigma=None,
        flatten=False,
    )

    assert light_curve.calls == [("remove_nans", {})]


def test_reduce_lightcurve_derives_pbjam_window_from_numax():
    """Numax should select the PBjam-compatible shorter flattening window."""
    light_curve = DummyLightCurve()

    mast.reduce_lightcurve(light_curve, numax=1000.0, exposure_time=100.0)

    assert light_curve.calls[-1] == ("flatten", {"window_length": 10001})


def test_reduce_lightcurve_infers_exposure_time_from_metadata():
    """Exposure time should be inferred from standard Lightkurve metadata."""
    light_curve = DummyLightCurve()

    mast.reduce_lightcurve(light_curve)

    expected = int(4e6 / 120.0)
    expected += int(expected % 2 == 0)
    assert light_curve.calls[-1] == ("flatten", {"window_length": expected})


def test_lightcurve_to_timeseries_converts_to_ppm():
    """Relative flux and uncertainties should be expressed in ppm together."""
    result = mast.lightcurve_to_timeseries(DummyLightCurve())

    assert np.allclose(result.time, [0.0, 1.0, 2.0])
    assert np.allclose(result.flux, [0.0, -1e5, 1e5])
    assert np.allclose(result.flux_err, [1e5, 1e5, 2e5])
    assert result.time_unit == "d"
    assert result.flux_unit == "ppm"


def test_lightcurve_to_timeseries_accepts_missing_uncertainties():
    """Light curves without flux errors should remain valid."""
    light_curve = DummyLightCurve()
    light_curve.flux_err = None

    result = mast.lightcurve_to_timeseries(light_curve, ppm=False)

    assert result.flux_err is None
    assert result.flux_unit is None


def test_load_lightcurve_orchestrates_separate_steps(monkeypatch):
    """The convenience loader should compose the public archive functions."""
    calls = []
    search = object()
    collection = object()
    reduced = DummyLightCurve()
    expected = object()

    monkeypatch.setattr(
        mast,
        "search_lightcurves",
        lambda target, **kwargs: calls.append(("search", target, kwargs)) or search,
    )
    monkeypatch.setattr(
        mast,
        "download_lightcurves",
        lambda value, **kwargs: calls.append(("download", value, kwargs))
        or collection,
    )
    monkeypatch.setattr(
        mast,
        "reduce_lightcurve",
        lambda value, **kwargs: calls.append(("reduce", value, kwargs)) or reduced,
    )
    monkeypatch.setattr(
        mast,
        "lightcurve_to_timeseries",
        lambda value, **kwargs: calls.append(("convert", value, kwargs)) or expected,
    )

    result = mast.load_lightcurve(
        "TIC 42",
        search_kwargs={"mission": "TESS"},
        download_dir="/cache",
        flatten=False,
        ppm=False,
    )

    assert result is expected
    assert calls[0] == ("search", "TIC 42", {"mission": "TESS"})
    assert calls[1] == ("download", search, {"download_dir": "/cache"})
    assert calls[2][0:2] == ("reduce", collection)
    assert calls[2][2]["flatten"] is False
    assert calls[3] == ("convert", reduced, {"ppm": False})


def test_load_lightcurve_reuses_search_exposure_time(monkeypatch):
    """The Lightkurve exptime constraint should also configure reduction."""
    captured = {}

    monkeypatch.setattr(mast, "search_lightcurves", lambda target, **kwargs: object())
    monkeypatch.setattr(
        mast,
        "download_lightcurves",
        lambda value, **kwargs: object(),
    )

    def reduce(value, **kwargs):
        """Record reduction options."""
        captured.update(kwargs)
        return DummyLightCurve()

    monkeypatch.setattr(mast, "reduce_lightcurve", reduce)
    monkeypatch.setattr(
        mast,
        "lightcurve_to_timeseries",
        lambda value, **kwargs: object(),
    )

    mast.load_lightcurve("KIC 123", search_kwargs={"exptime": 60})

    assert captured["exposure_time"] == 60
