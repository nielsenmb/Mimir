Power spectra
=============

Mimir computes regular-grid Lomb--Scargle periodograms with
`nifty-ls <https://github.com/flatironinstitute/nifty-ls>`_. The numerical
backend avoids the high-frequency bias that motivated moving away from the
standard Astropy fast implementation.

Basic use
---------

.. code-block:: python

   import numpy as np
   from mimir import power_spectrum

   time = np.arange(1000) * 120.0 / 86400.0
   flux = np.sin(2 * np.pi * 800e-6 * time * 86400.0)

   spectrum = power_spectrum(time, flux, oversampling=2)

   print(spectrum.frequency)
   print(spectrum.power_density)

An existing :class:`mimir.TimeSeries` is supplied through the explicit
``time_series`` argument. By default, input times are interpreted as days and
returned frequencies are in microhertz. Other Astropy-compatible time and
frequency units can be selected explicitly.

.. code-block:: python

   spectrum = power_spectrum(time_series=series, oversampling=2)

A target identifier can also be supplied directly. Mimir then obtains and
reduces the light curve through Lightkurve before calculating the spectrum:

.. code-block:: python

   spectrum = power_spectrum(
       "TIC 307210830",
       {"mission": "TESS", "author": "SPOC", "exptime": 120},
   )

Exactly one of the three input forms—array ``time`` with ``flux``,
``time_series``, or ``target``—must be selected. The same explicit convention
is supported by :func:`mimir.spectral_window`.

The target and flat ``mast_kwargs`` mapping may be positional, as shown above,
or named. Search filters go directly in ``mast_kwargs``. Optional download,
reduction, and conversion settings go in ``lightcurve_kwargs``, for example
``lightcurve_kwargs={"numax": 3500, "flatten": False}``. Legacy nested search
options are still accepted. Use :func:`mimir.search_lightcurves` to inspect
available products before downloading; see :doc:`mast` for the full workflow.

Threading
---------

By default, nifty-ls chooses how many threads to use. An explicit count can be
provided when its automatic choice oversubscribes the available physical CPU
cores:

.. code-block:: python

   spectrum = power_spectrum(time_series=series, nthreads=6)

When several targets are processed concurrently, use ``nthreads=1`` to avoid
nested parallelism, where every target starts its own pool of worker threads.

Frequency grid
--------------

The nominal Fourier spacing is ``1 / T``, where ``T`` is the duration of the
validated time series. ``oversampling`` divides that spacing, while
``nyquist_factor`` controls the upper limit relative to the Nyquist frequency
estimated from the median cadence. The zero-frequency bin is not returned.

For analyses that require every target on exactly the same grid, pass a
positive, increasing, regularly spaced array through ``frequency``. Its values
are interpreted in ``frequency_unit`` and are returned unchanged:

.. code-block:: python

   shared_frequency = np.arange(10.0, 5000.0, 0.1)  # microhertz
   spectrum = power_spectrum(
       time_series=series,
       frequency=shared_frequency,
       frequency_unit="uHz",
   )

An explicit grid cannot be combined with ``oversampling`` or
``nyquist_factor``. Mimir evaluates the requested bins with nifty-ls and uses a
fixed reference grid at positive multiples of ``1 / T`` through the
median-cadence Nyquist estimate, independent of the output spacing or band.
The ``oversampling`` metadata then
records the effective grid density, ``1 / (T * frequency_spacing)``, which may
differ between targets even though their returned frequency arrays are equal.

Normalization
-------------

Mimir fixes the density scale using a reference grid whose spacing is
``df_ref = 1 / T`` and whose frequencies are ``k * df_ref`` for
``k = 1, ..., floor(f_Nyquist / df_ref)``. Zero frequency is omitted because
the data are centred. The same reference is used regardless of output-grid
oversampling, spacing, starting frequency, or upper limit.

If ``R`` is the raw periodogram and ``V`` is the flux variance, the density is

.. math::

   S(f) = R(f)\,\frac{V}{\Delta f_{\rm ref}\sum_k R(k\Delta f_{\rm ref})}.

The default automatic grid (``oversampling=1, nyquist_factor=1``) is the
reference grid and its integrated density equals ``V`` by construction.
An oversampled grid estimates the integral numerically, a restricted band
represents only that band, and a coarse grid may miss narrow peaks. None of
these returned grids is rescaled to force its integral to equal ``V``.
``power`` is ``power_density * frequency_spacing`` using the **output**
spacing, not the reference spacing. Frequencies above the reference band do
not alter the normalization. With uncertainties, inverse-variance weights
define both the mean and variance.

This is a fixed variance-normalization convention. ``T`` is the last-minus-
first retained timestamp, so even for regular observations ``1/T`` is a
nominal Fourier spacing rather than exactly ``1/(N * cadence)``. For irregular
or gapped observations, the frequency basis is not orthogonal and the
median-cadence Nyquist frequency is a heuristic. The convention does not
establish exact Parseval reconstruction or independent output bins.

The default-spacing calculation reuses one backend evaluation. Oversampled
and explicit grids require an additional reference evaluation; ``nthreads``
applies to both. This reference is never made artificially coarse to match
a coarse requested output grid.

``amplitude`` is ``sqrt(2 * oversampling * power)`` and therefore has the same
units as the input flux. The oversampling factor compensates for narrower
frequency bins, keeping the semi-amplitude of a coherent sinusoid stable as the
grid is refined. Oversampled bins remain correlated and are not independent
Fourier measurements.

These definitions deliberately replace the dimensionally inconsistent legacy
amplitude calculation captured by the PBjam compatibility tests.

Spectral window
---------------

The sampling pattern can be inspected independently of the flux values:

.. code-block:: python

   from mimir import spectral_window

   window = spectral_window(time=time)
   print(window.effective_frequency_spacing)

The spectral window is the squared modulus of the discrete Fourier transform
of unit weights placed at the observation times. It is centred on zero,
normalized to unit height, and symmetric. Its integral estimates the effective
spacing between approximately independent frequency bins. For an uninterrupted
series this approaches ``1 / T``; gaps broaden the window and generally
increase the effective spacing.

By default, Mimir integrates over 100 nominal frequency spacings on either side
of zero and samples each nominal spacing ten times. ``half_width`` and
``oversampling`` can be changed when a more detailed window calculation is
needed.

API
---

.. autofunction:: mimir.power_spectrum

.. autoclass:: mimir.PowerSpectrum
   :members:

.. autofunction:: mimir.spectral_window

.. autofunction:: mimir.effective_frequency_spacing

.. autoclass:: mimir.SpectralWindow
   :members:
