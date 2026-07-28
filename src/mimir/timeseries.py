"""Time-series validation and preparation."""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray


class TimeSeries:
    """A validated astronomical time series.

    Parameters
    ----------
    time : array-like
        Sample times. Values are stored as a one-dimensional floating-point
        NumPy array.
    flux : array-like
        Flux measurements corresponding to ``time``.
    flux_err : array-like, optional
        One-sigma flux uncertainties. Finite uncertainties must be strictly
        positive.
    bad_mask : array-like of bool, optional
        Samples to remove, where ``True`` marks a bad sample.
    sort : bool, default=True
        Sort samples into ascending time order. If ``False``, unordered input
        raises :class:`ValueError`.
    invalid : {"drop", "raise"}, default="drop"
        Whether to remove or reject unmasked samples containing a non-finite
        time, flux, or flux uncertainty.
    time_unit : str, default="d"
        Unit label for the time values and derived cadence and duration.
    flux_unit : str, optional
        Unit label for the flux values and uncertainties.

    Notes
    -----
    Duplicate timestamps are rejected because Mimir cannot safely choose how
    independent measurements should be combined. Aggregate duplicates before
    constructing the object when that behaviour is appropriate.

    Duty cycle is estimated by comparing the retained sample count with the
    number expected across the time span at the median cadence. It is therefore
    an estimate for irregular data and is always bounded between zero and one.
    """

    def __init__(
        self,
        time: ArrayLike,
        flux: ArrayLike,
        flux_err: ArrayLike | None = None,
        *,
        bad_mask: ArrayLike | None = None,
        sort: bool = True,
        invalid: Literal["drop", "raise"] = "drop",
        time_unit: str = "d",
        flux_unit: str | None = None,
    ) -> None:
        """Validate and store the time-series samples."""
        if invalid not in {"drop", "raise"}:
            raise ValueError("invalid must be either 'drop' or 'raise'")

        time_values = _as_float_vector("time", time)
        flux_values = _as_float_vector("flux", flux)
        _require_matching_length(time_values, flux_values, "flux")

        error_values = None
        if flux_err is not None:
            error_values = _as_float_vector("flux_err", flux_err)
            _require_matching_length(time_values, error_values, "flux_err")

        input_size = time_values.size
        mask = _as_bad_mask(bad_mask, input_size)
        finite = np.isfinite(time_values) & np.isfinite(flux_values)
        if error_values is not None:
            finite &= np.isfinite(error_values)

        unmasked_invalid = ~finite & ~mask
        if invalid == "raise" and np.any(unmasked_invalid):
            count = int(np.count_nonzero(unmasked_invalid))
            raise ValueError(f"time series contains {count} non-finite sample(s)")

        keep = ~mask & finite
        time_values = time_values[keep]
        flux_values = flux_values[keep]
        if error_values is not None:
            error_values = error_values[keep]

        if time_values.size < 2:
            raise ValueError("time series must contain at least two valid samples")

        if error_values is not None and np.any(error_values <= 0.0):
            raise ValueError("flux_err values must be strictly positive")

        order = np.argsort(time_values, kind="stable")
        is_ordered = np.array_equal(order, np.arange(time_values.size))
        if not is_ordered and not sort:
            raise ValueError("time values must be strictly increasing")
        if not is_ordered:
            time_values = time_values[order]
            flux_values = flux_values[order]
            if error_values is not None:
                error_values = error_values[order]

        time_steps = np.diff(time_values)
        if np.any(time_steps == 0.0):
            raise ValueError("duplicate time values are not supported")

        self.time = time_values
        self.flux = flux_values
        self.flux_err = error_values
        self.time_unit = _validate_unit("time_unit", time_unit, optional=False)
        self.flux_unit = _validate_unit("flux_unit", flux_unit, optional=True)
        self.input_size = input_size
        self.n_removed = input_size - time_values.size

        self.cadence = float(np.median(time_steps))
        self.duration = float(time_values[-1] - time_values[0])
        expected_size = int(np.floor(self.duration / self.cadence + 0.5)) + 1
        self.duty_cycle = float(np.clip(time_values.size / expected_size, 0.0, 1.0))

    @property
    def n_samples(self) -> int:
        """Number of retained samples."""
        return int(self.time.size)

    @property
    def has_uncertainties(self) -> bool:
        """Whether flux uncertainties are available."""
        return self.flux_err is not None


def _as_float_vector(name: str, values: ArrayLike) -> NDArray[np.float64]:
    """Return an independent one-dimensional floating-point array."""
    try:
        array = np.array(values, dtype=float, copy=True)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must contain numerical values") from error
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def _require_matching_length(
    time: NDArray[np.float64],
    values: NDArray[np.float64],
    name: str,
) -> None:
    """Require an input vector to match the time vector."""
    if values.size != time.size:
        raise ValueError(
            f"{name} must have the same length as time "
            f"({values.size} != {time.size})"
        )


def _as_bad_mask(values: ArrayLike | None, size: int) -> NDArray[np.bool_]:
    """Validate a bad-sample mask."""
    if values is None:
        return np.zeros(size, dtype=bool)
    mask = np.asarray(values)
    if mask.ndim != 1:
        raise ValueError("bad_mask must be one-dimensional")
    if mask.size != size:
        raise ValueError(
            "bad_mask must have the same length as time "
            f"({mask.size} != {size})"
        )
    if not np.issubdtype(mask.dtype, np.bool_):
        raise TypeError("bad_mask must contain boolean values")
    return np.array(mask, dtype=bool, copy=True)


def _validate_unit(
    name: str,
    value: str | None,
    *,
    optional: bool,
) -> str | None:
    """Validate a unit label without imposing an array unit system."""
    if value is None and optional:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value
