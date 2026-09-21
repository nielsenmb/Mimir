"""Regression coverage for concise target input and product discovery."""

from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from mimir import (
    TimeSeries,
    as_timeseries,
    effective_frequency_spacing,
    mast,
    power_spectrum,
    search_lightcurves,
    spectral_window,
)


@pytest.mark.parametrize(
    "entrypoint",
    [power_spectrum, as_timeseries, spectral_window, effective_frequency_spacing],
)
def test_flat_target_workflow_preserves_filters_and_reduction(monkeypatch, entrypoint):
    """Resolve concise calls without altering archive or numerical options.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Loader replacement fixture.
    entrypoint : callable
        Public numerical entry point.
    """
    time = np.arange(256, dtype=float)
    series = TimeSeries(time, np.sin(time), time_unit="s")
    loader = Mock(return_value=series)
    monkeypatch.setattr(mast, "load_lightcurve", loader)
    filters = {"mission": "TESS", "author": "SPOC", "exptime": "short", "sector": 20}
    reduction = {"numax": 3500, "flatten": False}
    result = entrypoint("TIC 123", filters, lightcurve_kwargs=reduction)
    loader.assert_called_once_with("TIC 123", search_kwargs=filters, **reduction)
    assert filters == {
        "mission": "TESS", "author": "SPOC", "exptime": "short", "sector": 20,
    }
    assert reduction == {"numax": 3500, "flatten": False}
    if entrypoint is power_spectrum:
        direct = power_spectrum(time_series=series)
        np.testing.assert_allclose(result.power_density, direct.power_density)


def test_named_flat_and_legacy_loader_forms_match(monkeypatch):
    """Keep existing nested calls equivalent to the new named interface.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Loader replacement fixture.
    """
    loader = Mock(return_value=TimeSeries([0, 1], [1, 2]))
    monkeypatch.setattr(mast, "load_lightcurve", loader)
    as_timeseries(
        target="KIC 1", mast_kwargs={"mission": "Kepler"},
        lightcurve_kwargs={"numax": 3500},
    )
    as_timeseries(
        target="KIC 1",
        mast_kwargs={"search_kwargs": {"mission": "Kepler"}, "numax": 3500},
    )
    assert loader.call_args_list[0] == loader.call_args_list[1]
    as_timeseries("KIC 1")
    loader.assert_called_with("KIC 1", search_kwargs={})


@pytest.mark.parametrize(
    ("args", "kwargs"),
    [
        (("KIC 1", {}), {"mast_kwargs": {}}),
        (("KIC 1",), {"target": "KIC 2"}),
        (("KIC 1", []), {}),
        (("KIC 1", {}, {}), {}),
        (("KIC 1",), {"flux": [1, 2]}),
        (([0, 1],), {"time": [0, 1], "flux": [1, 2]}),
        (([0, 1], [1, 2]), {"lightcurve_kwargs": {}}),
        (("KIC 1",), {"mast_kwargs": {"mission": "TESS", "search_kwargs": {
            "mission": "Kepler"}}}),
        (("KIC 1",), {"mast_kwargs": {"numax": 3500}, "lightcurve_kwargs": {
            "numax": 3000}}),
        (("KIC 1",), {"lightcurve_kwargs": {"search_kwargs": {}}}),
    ],
)
def test_ambiguous_calls_fail_before_archive_access(monkeypatch, args, kwargs):
    """Reject duplicate inputs and malformed options before network access.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Loader replacement fixture.
    args : tuple
        Conflicting positional arguments.
    kwargs : dict
        Conflicting keyword arguments.
    """
    loader = Mock()
    monkeypatch.setattr(mast, "load_lightcurve", loader)
    with pytest.raises(TypeError):
        as_timeseries(*args, **kwargs)
    loader.assert_not_called()


def test_search_returns_product_metadata_without_downloading(monkeypatch):
    """Return the search result itself for inspection and native selection.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Lightkurve search replacement fixture.
    """
    products = [SimpleNamespace(author="SPOC"), SimpleNamespace(author="QLP")]
    search = Mock(return_value=products)
    monkeypatch.setattr(mast, "_import_lightkurve", lambda: SimpleNamespace(
        search_lightcurve=search,
    ))
    download = Mock()
    monkeypatch.setattr(mast, "download_lightcurves", download)
    result = search_lightcurves("TIC 1", {"mission": "TESS"}, exptime=120)
    assert result is products
    assert result[1].author == "QLP"
    search.assert_called_once_with("TIC 1", mission="TESS", exptime=120)
    download.assert_not_called()
    with pytest.raises(TypeError, match="duplicate search filters"):
        search_lightcurves("TIC 1", {"mission": "TESS"}, mission="Kepler")


@pytest.mark.parametrize("exptime", ["short", "long", "fast"])
def test_named_exposure_search_uses_observed_cadence(monkeypatch, exptime):
    """Keep named search cadence selectors out of numerical reduction settings.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Archive and conversion replacement fixture.
    exptime : str
        Valid Lightkurve search cadence selector.
    """
    search = Mock(return_value=object())
    monkeypatch.setattr(mast, "search_lightcurves", search)
    monkeypatch.setattr(mast, "download_lightcurves", Mock(return_value=object()))
    reduce = Mock(return_value=object())
    monkeypatch.setattr(mast, "reduce_lightcurve", reduce)
    monkeypatch.setattr(mast, "lightcurve_to_timeseries", Mock())
    mast.load_lightcurve("KIC 1", search_kwargs={"exptime": exptime})
    search.assert_called_once_with("KIC 1", exptime=exptime)
    assert reduce.call_args.kwargs["exposure_time"] is None
