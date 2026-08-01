"""Shared input resolution for Mimir's numerical entry points."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from mimir.timeseries import TimeSeries


def as_timeseries(
    time: ArrayLike | None = None,
    flux: ArrayLike | None = None,
    flux_err: ArrayLike | None = None,
    *,
    time_series: TimeSeries | None = None,
    target: str | None = None,
    mast_kwargs: Mapping[str, Any] | None = None,
    time_unit: str = "d",
    flux_unit: str | None = None,
) -> TimeSeries:
    """Return a validated time series from one explicit input form.

    Parameters
    ----------
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
        Options passed to :func:`mimir.load_lightcurve` when ``target`` is
        selected. Lightkurve search options belong in its nested
        ``search_kwargs`` mapping.
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

    Search and reduction options retain the loader structure:

    >>> series = as_timeseries(
    ...     target="KIC 8006161",
    ...     mast_kwargs={
    ...         "search_kwargs": {"mission": "Kepler", "exptime": 60},
    ...         "numax": 3500,
    ...     },
    ... )
    """
    return _resolve_timeseries_input(
        time=time,
        flux=flux,
        flux_err=flux_err,
        time_series=time_series,
        target=target,
        mast_kwargs=mast_kwargs,
        time_unit=time_unit,
        flux_unit=flux_unit,
        allow_time_only=False,
    )


def _resolve_timeseries_input(
    time: ArrayLike | None,
    flux: ArrayLike | None,
    flux_err: ArrayLike | None,
    *,
    time_series: TimeSeries | None,
    target: str | None,
    mast_kwargs: Mapping[str, Any] | None,
    time_unit: str,
    flux_unit: str | None,
    allow_time_only: bool,
) -> TimeSeries:
    """Resolve mutually exclusive array, object, and target inputs."""
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
        return _load_mast_target(target, mast_kwargs)

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


def _load_mast_target(
    target: str,
    mast_kwargs: Mapping[str, Any] | None,
) -> TimeSeries:
    """Load one target lazily so Lightkurve remains an optional dependency."""
    from mimir.mast import load_lightcurve

    kwargs = {} if mast_kwargs is None else dict(mast_kwargs)
    return load_lightcurve(target, **kwargs)
