"""Shared input resolution for Mimir's numerical entry points."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from mimir.timeseries import TimeSeries


def as_timeseries(
    *args: Any,
    time: ArrayLike | None = None,
    flux: ArrayLike | None = None,
    flux_err: ArrayLike | None = None,
    time_series: TimeSeries | None = None,
    target: str | None = None,
    mast_kwargs: Mapping[str, Any] | None = None,
    lightcurve_kwargs: Mapping[str, Any] | None = None,
    time_unit: str = "d",
    flux_unit: str | None = None,
) -> TimeSeries:
    """Return a validated time series from one explicit input form.

    Parameters
    ----------
    *args
        Target name and optional flat search-filter mapping, for example
        ``("KIC 8006161", {"mission": "Kepler", "exptime": 60})``.
        Positional sample arrays remain supported; do not repeat named inputs.
    time : array-like, optional
        Sample times. ``flux`` is required when this input form is selected.
    flux : array-like, optional
        Flux measurements corresponding to ``time``.
    flux_err : array-like, optional
        Positive one-sigma flux uncertainties corresponding to ``time``.
    time_series : TimeSeries, optional
        Existing validated time series. Values of other types are rejected.
    target : str, optional
        Target name or identifier understood by Lightkurve. A target triggers
        MAST retrieval through :func:`mimir.load_lightcurve`.
    mast_kwargs : mapping, optional
        Flat Lightkurve search filters, such as ``mission``, ``author``,
        ``exptime``, ``quarter``, or ``sector``. No nesting is required.
        The legacy nested ``search_kwargs`` form remains supported.
    lightcurve_kwargs : mapping, optional
        Download, reduction, and conversion options for
        :func:`mimir.load_lightcurve`, such as ``download_dir``, ``numax``,
        ``flatten``, or ``ppm``. Only valid with target input.
    time_unit : str, default="d"
        Unit label for array-like sample times.
    flux_unit : str, optional
        Unit label for array-like flux values.

    Returns
    -------
    TimeSeries
        Downloaded, reused, or newly validated time series.

    Raises
    ------
    TypeError
        If no input form is selected, more than one input form is selected,
        or an input has the wrong type.

    Examples
    --------
    Reuse an existing validated object explicitly:

    >>> series = as_timeseries(time_series=validated_series)

    Download a target that Lightkurve can resolve:

    >>> series = as_timeseries(target="KIC 8006161")

    Use a target and a flat mapping for a quick download:

    >>> series = as_timeseries("KIC 8006161", {"mission": "Kepler", "exptime": 60})

    Configure reduction separately when needed:

    >>> series = as_timeseries(
    ...     target="KIC 8006161",
    ...     mast_kwargs={"mission": "Kepler", "exptime": 60},
    ...     lightcurve_kwargs={"numax": 3500},
    ... )
    """
    return _resolve_timeseries_input(
        *args,
        time=time,
        flux=flux,
        flux_err=flux_err,
        time_series=time_series,
        target=target,
        mast_kwargs=mast_kwargs,
        lightcurve_kwargs=lightcurve_kwargs,
        time_unit=time_unit,
        flux_unit=flux_unit,
        allow_time_only=False,
    )


def _resolve_timeseries_input(
    *args: Any,
    time: ArrayLike | None,
    flux: ArrayLike | None,
    flux_err: ArrayLike | None,
    time_series: TimeSeries | None,
    target: str | None,
    mast_kwargs: Mapping[str, Any] | None,
    lightcurve_kwargs: Mapping[str, Any] | None,
    time_unit: str,
    flux_unit: str | None,
    allow_time_only: bool,
) -> TimeSeries:
    """Resolve positional convenience calls and distinct named inputs.

    Parameters
    ----------
    *args
        Target and optional filters, or positional sample arrays.
    time, flux, flux_err : array-like or None
        Explicit sample arrays.
    time_series : TimeSeries or None
        Existing validated object.
    target : str or None
        Explicit target name.
    mast_kwargs : mapping or None
        Search filters or legacy loader options.
    lightcurve_kwargs : mapping or None
        Download, reduction, and conversion settings.
    time_unit, flux_unit : str or None
        Unit labels for array input.
    allow_time_only : bool
        Whether a sampling-window call may omit flux.

    Returns
    -------
    TimeSeries
        Reused, downloaded, or newly validated measurements.
    """
    if args:
        if isinstance(args[0], str):
            if len(args) > 2:
                raise TypeError("target input accepts a name and optional mast_kwargs")
            if target is not None or time is not None or time_series is not None:
                raise TypeError("select exactly one input form; target was positional")
            target = args[0]
            if len(args) == 2:
                if mast_kwargs is not None:
                    raise TypeError("mast_kwargs was supplied twice")
                mast_kwargs = _options_mapping("mast_kwargs", args[1])
        else:
            names = ("time",) if allow_time_only else ("time", "flux", "flux_err")
            if len(args) > len(names):
                raise TypeError(f"array input accepts at most {len(names)} arguments")
            values = {"time": time, "flux": flux, "flux_err": flux_err}
            for name, value in zip(names, args, strict=False):
                if values[name] is not None:
                    raise TypeError(f"{name} was supplied twice")
                values[name] = value
            time, flux, flux_err = values.values()
    if lightcurve_kwargs is not None and target is None:
        raise TypeError("lightcurve_kwargs can only be used with target")
    selected = sum(value is not None for value in (time, time_series, target))
    if selected != 1:
        raise TypeError(
            "select exactly one input form: time with flux, time_series, or target"
        )

    if time_series is not None:
        if not isinstance(time_series, TimeSeries):
            raise TypeError("time_series must be a TimeSeries object")
        _reject_arguments(flux, flux_err, mast_kwargs, "time_series")
        return time_series

    if target is not None:
        if not isinstance(target, str):
            raise TypeError("target must be a string")
        if not target.strip():
            raise ValueError("target must not be empty")
        _reject_arguments(flux, flux_err, None, "target")
        if lightcurve_kwargs is None:
            return _load_mast_target(target, mast_kwargs)
        return _load_mast_target(target, mast_kwargs, lightcurve_kwargs)

    if mast_kwargs is not None:
        raise TypeError("mast_kwargs can only be used with target")
    if flux is None:
        if not allow_time_only:
            raise TypeError("flux is required when time is supplied")
        values = _one_dimensional_time(time)
        flux = np.zeros(values.size, dtype=float)
        time = values
    return TimeSeries(
        time=time,
        flux=flux,
        flux_err=flux_err,
        time_unit=time_unit,
        flux_unit=flux_unit,
    )


def _one_dimensional_time(time: ArrayLike | None) -> np.ndarray:
    """Return numerical one-dimensional sample times."""
    try:
        values = np.asarray(time, dtype=float)
    except (TypeError, ValueError) as error:
        raise TypeError("time must contain numerical values") from error
    if values.ndim != 1:
        raise ValueError("time must be one-dimensional")
    return values


def _reject_arguments(
    flux: ArrayLike | None,
    flux_err: ArrayLike | None,
    mast_kwargs: Mapping[str, Any] | None,
    input_name: str,
) -> None:
    """Reject arguments that do not apply to the selected input form."""
    if flux is not None or flux_err is not None:
        raise TypeError(f"flux and flux_err must be omitted with {input_name}")
    if mast_kwargs is not None:
        raise TypeError(f"mast_kwargs must be omitted with {input_name}")


def _options_mapping(name: str, value: Mapping[str, Any] | None) -> dict[str, Any]:
    """Copy an optional mapping of keyword options.

    Parameters
    ----------
    name : str
        Parameter name to use in validation errors.
    value : mapping or None
        Options with string keys.

    Returns
    -------
    dict
        Independent top-level copy of the options.
    """
    if value is None:
        return {}
    if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
        raise TypeError(f"{name} must be a mapping with string keys")
    return dict(value)


def _load_mast_target(
    target: str,
    mast_kwargs: Mapping[str, Any] | None,
    lightcurve_kwargs: Mapping[str, Any] | None = None,
) -> TimeSeries:
    """Load a target using flat filters and separate light-curve options.

    Parameters
    ----------
    target : str
        Target name understood by Lightkurve.
    mast_kwargs : mapping or None
        Flat search filters, or the legacy mapping of loader options.
    lightcurve_kwargs : mapping or None, optional
        Download, reduction, and conversion settings.

    Returns
    -------
    TimeSeries
        Downloaded and reduced time series.

    Raises
    ------
    TypeError
        If options are invalid or supplied through more than one mapping.
    """
    from mimir.mast import load_lightcurve

    filters = _options_mapping("mast_kwargs", mast_kwargs)
    options = _options_mapping("lightcurve_kwargs", lightcurve_kwargs)
    if "target" in options or "search_kwargs" in options:
        raise TypeError("put the target in target and search filters in mast_kwargs")
    nested = _options_mapping("search_kwargs", filters.pop("search_kwargs", None))
    legacy_names = {
        "download_dir", "outlier_sigma", "flatten", "flatten_window_length",
        "exposure_time", "numax", "normalize", "ppm",
    }
    for name in legacy_names & filters.keys():
        if name in options:
            raise TypeError(f"{name} was supplied in both option mappings")
        options[name] = filters.pop(name)
    if duplicate := nested.keys() & filters.keys():
        raise TypeError(f"duplicate search filters: {', '.join(sorted(duplicate))}")
    if "target" in filters:
        raise TypeError("put the target in target, not mast_kwargs")
    filters.update(nested)
    return load_lightcurve(target, search_kwargs=filters, **options)
