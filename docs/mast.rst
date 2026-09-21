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

Compute a spectrum using a target and a flat dictionary of search filters:

.. code-block:: python

   from mimir import power_spectrum

   target = "KIC 8006161"
   mast_kwargs = {"mission": "Kepler", "author": "Kepler", "exptime": 60}
   spectrum = power_spectrum(target, mast_kwargs)

If Lightkurve can resolve the name using its defaults, the only required
argument is the target:

.. code-block:: python

   spectrum = power_spectrum("KIC 8006161")

``mast_kwargs`` contains Lightkurve search filters: ``mission``, ``author``,
``exptime``, ``quarter``, ``sector``, and other supported search keywords.
There is no nested ``search_kwargs`` dictionary to construct. Named arguments
work too: ``power_spectrum(target=target, mast_kwargs=mast_kwargs)``.

The convenience call downloads all matching products. With multiple pipelines
or cadences available, inspect the list first and constrain the search to the
products you intend to combine.

Inspect and select products
---------------------------

Search without downloading any light curves:

.. code-block:: python

   from mimir import search_lightcurves

   products = search_lightcurves(target)
   print(products)
   products.table  # Full product metadata; displays as a table in a notebook

This returns Lightkurve's native ``SearchResult``. It includes supported
high-level science light-curve products, not an inventory of every data type
stored at MAST. ``author`` selects the data-producing pipeline; omitting it
allows all authors supported by the search. The same flat filters can be used
for discovery and calculation:

.. code-block:: python

   products = search_lightcurves(target, mast_kwargs)
   # Equivalent keyword form:
   products = search_lightcurves(target, **mast_kwargs)

Inspect the rows and select the desired subset before downloading. For example,
after deciding that the first row is the product you want:

.. code-block:: python

   from mimir import download_lightcurves, reduce_lightcurve, lightcurve_to_timeseries

   selected = products[:1]  # Choose indices using the displayed metadata
   collection = download_lightcurves(selected)
   lightcurve = reduce_lightcurve(collection, numax=3500)
   series = lightcurve_to_timeseries(lightcurve)
   spectrum = power_spectrum(time_series=series)

Further details on search filters and result selection are in the
`Lightkurve search documentation <https://lightkurve.github.io/lightkurve/reference/api/lightkurve.search_lightcurve.html>`_.

Download and reduction options
------------------------------

For optional processing settings, use ``lightcurve_kwargs``:

.. code-block:: python

   spectrum = power_spectrum(
       target,
       mast_kwargs,
       lightcurve_kwargs={"numax": 3500, "outlier_sigma": 5.0},
   )

These are the download, reduction, and conversion settings accepted by
:func:`mimir.load_lightcurve`, such as ``download_dir``, ``flatten``, and
``ppm``. Lightkurve's named ``exptime`` selectors (``"short"``, ``"long"``, and
``"fast"``) are search filters; Mimir infers the numerical cadence from the
downloaded light curve when one of those selectors is used.

Existing nested ``mast_kwargs={"search_kwargs": {...}, "numax": ...}`` calls
remain supported. New code can use the flat form. Supplying the same option
through two mappings raises an error instead of silently overriding it.

The :func:`mimir.as_timeseries` wrapper provides the same explicit input handling
without computing a spectrum. Explicit :func:`mimir.load_lightcurve` calls
remain available when archive retrieval and numerical calculation should be
separate.

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

Time conversion and masked measurements
---------------------------------------

:func:`mimir.lightcurve_to_timeseries` uses Astropy time arithmetic rather
than the displayed ``Time.value``. Absolute ``Time`` columns are converted to
elapsed times since JD 2451545.0 in their input time scale; their output no
longer retains the original BKJD, BTJD, Unix, or JD display offset. Thus
changing only the display format cannot change cadence, frequency, or PSD
units. ``time_unit="s"`` performs an actual conversion to seconds. The
origin is JD 2451545.0 in the input scale, not a conversion of all inputs to
TT or TDB. UTC differences include leap seconds through Astropy arithmetic.

``TimeDelta`` and time ``Quantity`` columns are converted directly without
changing their origin. Unitless numerical columns are assumed to be days.
The exposure-time fallback uses the same conversion to seconds before
estimating cadence. Explicit exposure metadata still takes precedence.

Masks on time, flux, and uncertainty columns are applied before computing
the ppm normalization median. Masked values cannot re-enter as finite flux
outliers. The input light curve is not modified, and the returned time
series retains its input and removal counts.

API
---

.. autofunction:: mimir.search_lightcurves

.. autofunction:: mimir.download_lightcurves

.. autofunction:: mimir.reduce_lightcurve

.. autofunction:: mimir.lightcurve_to_timeseries

.. autofunction:: mimir.load_lightcurve

.. autofunction:: mimir.as_timeseries
