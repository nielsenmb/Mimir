"""Optional MAST access and Lightkurve-based light-curve reduction.

Archive access is kept separate from the numerical spectrum implementation so
that Lightkurve remains an optional dependency.
"""

from __future__ import annotations

from importlib import import_module
from operator import index
from os import PathLike
from typing import Any

import numpy as np

from mimir.timeseries import TimeSeries


def search_lightcurves(target: str, **search_kwargs: Any) -> Any:
    """Search MAST for light curves using Lightkurve.

    Parameters
    ----------
    target : str
        Target identifier understood by
        :func:`lightkurve.search_lightcurve`.
    **search_kwargs
        Additional search constraints, such as ``mission``, ``author``, or
        ``exptime``.

    Returns
    -------
    lightkurve.SearchResult
        Non-empty Lightkurve search result.

    Raises
    ------
    ImportError
        If the optional Lightkurve dependency is not installed.
    LookupError
        If the search returns no matching light curves.
    ValueError
        If ``target`` is not a non-empty string.
    """
    if not isinstance(target, str) or not target.strip():
        raise ValueError("target must be a non-empty string")

    lightkurve = _import_lightkurve()
    result = lightkurve.search_lightcurve(target, **search_kwargs)
    if _is_empty(result):
        raise LookupError(f"no light curves found for {target!r}")
    return result


def download_lightcurves(
    search_result: Any,
    *,
    download_dir: str | PathLike[str] | None = None,
) -> Any:
    """Download every light curve in a Lightkurve search result.

    Parameters
    ----------
    search_result : lightkurve.SearchResult
        Search result returned by :func:`search_lightcurves`.
    download_dir : path-like, optional
        Directory used by Lightkurve as its download cache.

    Returns
    -------
    lightkurve.LightCurveCollection
        Non-empty collection of downloaded light curves.

    Raises
    ------
    LookupError
        If the search result or downloaded collection is empty.
    TypeError
        If ``search_result`` does not provide ``download_all``.
    """
    if _is_empty(search_result):
        raise LookupError("cannot download an empty search result")
    if not hasattr(search_result, "download_all"):
        raise TypeError("search_result must provide a download_all method")

    collection = search_result.download_all(download_dir=download_dir)
    if _is_empty(collection):
        raise LookupError("Lightkurve did not download any light curves")
    return collection


def reduce_lightcurve(
    light_curve: Any,
    *,
    outlier_sigma: float | None = 5.0,
    flatten: bool = True,
    flatten_window_length: int | None = None,
    exposure_time: float | None = None,
    numax: float | None = None,
    normalize: bool = True,
) -> Any:
    """Stitch and perform basic Lightkurve light-curve reduction.

    Parameters
    ----------
    light_curve : lightkurve.LightCurve or lightkurve.LightCurveCollection
        A light curve or collection of light curves. Collections are stitched
        before reduction.
    outlier_sigma : float or None, default=5.0
        Sigma threshold passed to ``remove_outliers``. Set to ``None`` to skip
        outlier removal.
    flatten : bool, default=True
        Whether to remove long-period trends with ``flatten``.
    flatten_window_length : int, optional
        Savitzky--Golay window length in cadences. Even values are advanced to
        the next odd integer. If omitted, a PBjam-compatible window is derived
        from ``exposure_time`` and optionally ``numax``.
    exposure_time : float, optional
        Exposure time in seconds used to derive the flattening window. Mimir
        attempts to infer it from the Lightkurve object if omitted.
    numax : float, optional
        Frequency of maximum oscillation power in microhertz. When provided,
        it selects the shorter PBjam-compatible flattening window.
    normalize : bool, default=True
        Whether to call Lightkurve's ``normalize`` before cleaning.

    Returns
    -------
    lightkurve.LightCurve
        Reduced light curve.

    Raises
    ------
    TypeError
        If the input does not provide the required Lightkurve methods.
    ValueError
        If a numerical reduction setting is invalid or the exposure time
        cannot be determined when flattening is requested.
    """
    reduced = _stitch(light_curve)

    if normalize:
        reduced = _call_method(reduced, "normalize")
    reduced = _call_method(reduced, "remove_nans")

    if outlier_sigma is not None:
        sigma = _positive_finite("outlier_sigma", outlier_sigma)
        reduced = _call_method(reduced, "remove_outliers", sigma=sigma)

    if flatten:
        if flatten_window_length is None:
            cadence = exposure_time
            if cadence is None:
                cadence = _infer_exposure_time(reduced)
            flatten_window_length = _pbjam_window_length(
                cadence,
                numax=numax,
            )
        else:
            flatten_window_length = _odd_window_length(flatten_window_length)
        reduced = _call_method(
            reduced,
            "flatten",
            window_length=flatten_window_length,
        )

    return reduced


def lightcurve_to_timeseries(
    light_curve: Any,
    *,
    ppm: bool = True,
    time_unit: str = "d",
) -> TimeSeries:
    """Convert a Lightkurve light curve into a Mimir time series.

    Parameters
    ----------
    light_curve : lightkurve.LightCurve
        Light curve containing ``time`` and ``flux`` columns and, optionally,
        ``flux_err``.
    ppm : bool, default=True
        Convert flux to relative parts per million using its finite median.
        Flux uncertainties are scaled by the same factor.
    time_unit : str, default="d"
        Unit label for the extracted time values.

    Returns
    -------
    TimeSeries
        Validated NumPy-based Mimir time series.

    Raises
    ------
    TypeError
        If the object does not provide time and flux columns.
    ValueError
        If ppm conversion is requested for a zero or non-finite median flux.
    """
    if not hasattr(light_curve, "time") or not hasattr(light_curve, "flux"):
        raise TypeError("light_curve must provide time and flux columns")

    time = _column_values(light_curve.time)
    flux = _column_values(light_curve.flux)
    flux_err_column = getattr(light_curve, "flux_err", None)
    flux_err = (
        None if flux_err_column is None else _column_values(flux_err_column)
    )

    flux_unit = None
    if ppm:
        median = float(np.nanmedian(flux))
        if not np.isfinite(median) or median == 0.0:
            raise ValueError("cannot convert a zero or non-finite median flux to ppm")
        scale = 1e6 / abs(median)
        flux = (flux / median - 1.0) * 1e6
        if flux_err is not None:
            flux_err = flux_err * scale
        flux_unit = "ppm"

    return TimeSeries(
        time,
        flux,
        flux_err,
        time_unit=time_unit,
        flux_unit=flux_unit,
    )


def load_lightcurve(
    target: str,
    *,
    search_kwargs: dict[str, Any] | None = None,
    download_dir: str | PathLike[str] | None = None,
    outlier_sigma: float | None = 5.0,
    flatten: bool = True,
    flatten_window_length: int | None = None,
    exposure_time: float | None = None,
    numax: float | None = None,
    normalize: bool = True,
    ppm: bool = True,
) -> TimeSeries:
    """Search, download, reduce, and convert a MAST light curve.

    Parameters
    ----------
    target : str
        Target identifier understood by Lightkurve.
    search_kwargs : dict, optional
        Constraints passed to :func:`search_lightcurves`.
    download_dir : path-like, optional
        Directory used by Lightkurve as its download cache.
    outlier_sigma : float or None, default=5.0
        Sigma threshold used for outlier removal.
    flatten : bool, default=True
        Whether to remove long-period trends.
    flatten_window_length : int, optional
        Explicit Lightkurve flattening window in cadences.
    exposure_time : float, optional
        Exposure time in seconds used to derive a flattening window.
    numax : float, optional
        Frequency of maximum oscillation power in microhertz.
    normalize : bool, default=True
        Whether to normalize the stitched light curve.
    ppm : bool, default=True
        Whether to return relative flux in parts per million.

    Returns
    -------
    TimeSeries
        Downloaded and validated Mimir time series.
    """
    kwargs = {} if search_kwargs is None else dict(search_kwargs)
    if exposure_time is None and "exptime" in kwargs:
        exposure_time = kwargs["exptime"]
    result = search_lightcurves(target, **kwargs)
    collection = download_lightcurves(result, download_dir=download_dir)
    reduced = reduce_lightcurve(
        collection,
        outlier_sigma=outlier_sigma,
        flatten=flatten,
        flatten_window_length=flatten_window_length,
        exposure_time=exposure_time,
        numax=numax,
        normalize=normalize,
    )
    return lightcurve_to_timeseries(reduced, ppm=ppm)


def _import_lightkurve() -> Any:
    """Import and return the optional Lightkurve dependency."""
    try:
        return import_module("lightkurve")
    except ImportError as error:
        raise ImportError(
            "MAST access requires Lightkurve; install mimir-astro[mast]"
        ) from error


def _is_empty(value: Any) -> bool:
    """Return whether an archive result is absent or has zero length."""
    if value is None:
        return True
    try:
        return len(value) == 0
    except TypeError:
        return False


def _stitch(light_curve: Any) -> Any:
    """Return one light curve, stitching a collection when necessary."""
    if hasattr(light_curve, "time") and hasattr(light_curve, "flux"):
        return light_curve
    if not hasattr(light_curve, "stitch"):
        raise TypeError("light_curve must be a LightCurve or LightCurveCollection")
    stitched = light_curve.stitch()
    if stitched is None:
        raise LookupError("Lightkurve could not stitch the downloaded light curves")
    return stitched


def _call_method(value: Any, name: str, **kwargs: Any) -> Any:
    """Call a required Lightkurve reduction method and return its result."""
    method = getattr(value, name, None)
    if method is None:
        raise TypeError(f"light_curve must provide a {name} method")
    result = method(**kwargs)
    if result is None:
        raise RuntimeError(f"Lightkurve {name} returned None")
    return result


def _positive_finite(name: str, value: float) -> float:
    """Validate and return a positive finite number."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a positive finite number") from error
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return result


def _odd_window_length(value: int) -> int:
    """Validate a window length and advance even values to the next odd value."""
    if isinstance(value, bool):
        raise TypeError("flatten_window_length must be a positive integer")
    try:
        result = index(value)
    except TypeError as error:
        raise TypeError("flatten_window_length must be a positive integer") from error
    if result < 1:
        raise ValueError("flatten_window_length must be a positive integer")
    return result if result % 2 else result + 1


def _pbjam_window_length(
    exposure_time: float,
    *,
    numax: float | None,
) -> int:
    """Return the PBjam-compatible flattening window in cadences."""
    cadence = _positive_finite("exposure_time", exposure_time)
    if numax is None:
        window_seconds = 4e6
    else:
        numax_value = _positive_finite("numax", numax)
        window_seconds = 1e9 / numax_value
    return _odd_window_length(max(1, int(window_seconds / cadence)))


def _infer_exposure_time(light_curve: Any) -> float:
    """Infer exposure time in seconds from metadata or time samples."""
    meta = getattr(light_curve, "meta", {}) or {}
    if "EXPTIME" in meta:
        return _positive_finite("EXPTIME", meta["EXPTIME"])
    if "TIMEDEL" in meta:
        return _positive_finite("TIMEDEL", meta["TIMEDEL"]) * 86400.0

    if not hasattr(light_curve, "time"):
        raise ValueError(
            "exposure_time is required when it cannot be inferred from the light curve"
        )
    time = _column_values(light_curve.time)
    finite = np.sort(time[np.isfinite(time)])
    if finite.size < 2:
        raise ValueError(
            "exposure_time is required when it cannot be inferred from the light curve"
        )
    return _positive_finite("inferred exposure_time", np.median(np.diff(finite))) * (
        86400.0
    )


def _column_values(column: Any) -> np.ndarray:
    """Return floating-point NumPy values from a Lightkurve column."""
    values = getattr(column, "value", column)
    try:
        return np.asarray(values, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("Lightkurve columns must contain numerical values") from error
