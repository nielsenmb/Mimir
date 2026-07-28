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
from mimir.window import (
    SpectralWindow,
    effective_frequency_spacing,
    spectral_window,
)

try:
    __version__ = version("mimir-astro")
except PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = [
    "PowerSpectrum",
    "SpectralWindow",
    "TimeSeries",
    "__version__",
    "download_lightcurves",
    "effective_frequency_spacing",
    "lightcurve_to_timeseries",
    "load_lightcurve",
    "power_spectrum",
    "reduce_lightcurve",
    "search_lightcurves",
    "spectral_window",
]

