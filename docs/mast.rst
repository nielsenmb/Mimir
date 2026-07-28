MAST light curves
=================

Mimir keeps archive access separate from time-series validation and spectrum
calculation. Lightkurve handles MAST product discovery, download caching, and
mission-specific light-curve objects; Mimir provides a small, explicit layer
around those operations.

Install the optional dependency with:

.. code-block:: console

   python -m pip install "mimir-astro[mast]"

Basic use
---------

.. code-block:: python

   from mimir import load_lightcurve, power_spectrum

   series = load_lightcurve(
       "KIC 8006161",
       search_kwargs={"mission": "Kepler", "exptime": 60},
       numax=3500,
   )
   spectrum = power_spectrum(series)

The convenience loader composes four operations that are also available
independently:

#. :func:`mimir.search_lightcurves` searches MAST.
#. :func:`mimir.download_lightcurves` downloads all matches.
#. :func:`mimir.reduce_lightcurve` stitches and performs basic reduction.
#. :func:`mimir.lightcurve_to_timeseries` converts to validated NumPy arrays.

Reduction
---------

The default reduction follows PBjam's established sequence: normalize, remove
non-finite values, remove outliers, and flatten. Each optional operation can be
configured or disabled. The flattening window can be supplied directly; when
it is omitted, Mimir uses the existing PBjam heuristic based on exposure time
and, when available, ``numax``.

Network access is mocked in the ordinary test suite. Mimir's CI therefore does
not depend on MAST availability or create repeated archive requests.

API
---

.. autofunction:: mimir.search_lightcurves

.. autofunction:: mimir.download_lightcurves

.. autofunction:: mimir.reduce_lightcurve

.. autofunction:: mimir.lightcurve_to_timeseries

.. autofunction:: mimir.load_lightcurve
