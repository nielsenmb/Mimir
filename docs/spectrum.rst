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
       target="TIC 307210830",
       mast_kwargs={"search_kwargs": {"mission": "TESS", "exptime": 120}},
   )

Exactly one of the three input forms—array ``time`` with ``flux``,
``time_series``, or ``target``—must be selected. The same explicit convention
is supported by :func:`mimir.spectral_window`.

Frequency grid
--------------

The nominal Fourier spacing is ``1 / T``, where ``T`` is the duration of the
validated time series. ``oversampling`` divides that spacing, while
``nyquist_factor`` controls the upper limit relative to the Nyquist frequency
estimated from the median cadence. The zero-frequency bin is not returned.

Normalization
-------------

``power`` is a one-sided per-bin quantity. At the default
``nyquist_factor=1``, its sum equals the variance of the input flux.
``power_density`` is ``power / frequency_spacing``, so integrating it through
Nyquist gives the same result. Truncating the returned band retains only the
power represented there. Extending above Nyquist returns aliases without using
them to renormalize, and therefore dilute, the physical one-sided spectrum.
When flux uncertainties are provided, inverse-variance weights define both the
mean and the variance.

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

