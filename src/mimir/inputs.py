"""Shared input coercion for Mimir's numerical entry points."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from numpy.typing import ArrayLike

from mimir.timeseries import TimeSeries


def as_timeseries(
    source: TimeSeries | str | ArrayLike,
    flux: ArrayLike | None = None,
    flux_err: ArrayLike | None = None,
    *,
    mast_kwargs: Mapping[str, Any] | None = None,
    time_unit: str = "d",
    flux_unit: str | None = None,
) -> TimeSeries:
    """Return a validated time series from arrays, an object, or a MAST target.

    Parameters
    ----------
    source : TimeSeries, str, or array-like
        Existing time series, a target identifier understood by Lightkurve, or
        sample times. A string triggers MAST retrieval through
        :func:`mimir.load_lightcurve`.
    flux : array-like, optional
        Flux measurements. Required for array-like ``source`` and omitted for
        a ``TimeSeries`` or target identifier.
    flux_err : array-like, optional
        Positive one-sigma flux uncertainties. Omitted for a ``TimeSeries`` or
        target identifier.
    mast_kwargs : mapping, optional
        Options passed to :func:`mimir.load_lightcurve` when ``source`` is a
        target identifier. Lightkurve search options belong in its nested
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
        If arguments from different input forms are mixed, or array-like input
        is supplied without flux values.

    Examples
    --------
    A target name alone is sufficient when Lightkurve can resolve it:

    >>> series = as_timeseries("KIC 8006161")

    Search and reduction options retain the existing loader structure:

    >>> series = as_timeseries(
    ...     "KIC 8006161",
    ...     mast_kwargs={
    ...         "search_kwargs": {"mission": "Kepler", "exptime": 60},
    ...         "numax": 3500,
    ...     },
    ... )
    """
    if isinstance(source, TimeSeries):
        _reject_mixed_input(flux, flux_err, mast_kwargs, "a TimeSeries")
        return source

    if isinstance(source, str):
        _reject_mixed_input(flux, flux_err, None, "a target identifier")
        return _load_mast_target(source, mast_kwargs)

    if mast_kwargs is not None:
        raise TypeError("mast_kwargs can only be used with a target identifier")
    if flux is None:
        raise TypeError("flux is required when source contains sample times")
    return TimeSeries(
        time=source,
        flux=flux,
        flux_err=flux_err,
        time_unit=time_unit,
        flux_unit=flux_unit,
    )


def _reject_mixed_input(
    flux: ArrayLike | None,
    flux_err: ArrayLike | None,
    mast_kwargs: Mapping[str, Any] | None,
    source_description: str,
) -> None:
    """Reject arguments that do not apply to the selected input form."""
    if flux is not None or flux_err is not None:
        raise TypeError(
            f"flux and flux_err must be omitted when source is {source_description}"
        )
    if mast_kwargs is not None:
        raise TypeError(
            f"mast_kwargs must be omitted when source is {source_description}"
        )


def _load_mast_target(
    target: str,
    mast_kwargs: Mapping[str, Any] | None,
) -> TimeSeries:
    """Load one target lazily so Lightkurve remains an optional dependency."""
    from mimir.mast import load_lightcurve

    kwargs = {} if mast_kwargs is None else dict(mast_kwargs)
    return load_lightcurve(target, **kwargs)
