"""Shared time-series and power-spectrum tools for asteroseismology."""

from importlib.metadata import PackageNotFoundError, version

from mimir.mast import (
    download_lightcurves,
    lightcurve_to_timeseries,
    load_lightcurve,
    reduce_lightcurve,
    search_lightcurves,
)
from mimir.spectrum import PowerSpectrum, power_spectrum
from mimir.timeseries import TimeSeries

try:
    __version__ = version("mimir-astro")
except PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = [
    "PowerSpectrum",
    "TimeSeries",
    "__version__",
    "download_lightcurves",
    "lightcurve_to_timeseries",
    "load_lightcurve",
    "power_spectrum",
    "reduce_lightcurve",
    "search_lightcurves",
]
