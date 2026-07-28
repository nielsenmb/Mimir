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

An existing :class:`mimir.TimeSeries` can be supplied instead of separate
arrays. By default, input times are interpreted as days and returned
frequencies are in microhertz. Other Astropy-compatible time and frequency
units can be selected explicitly.

Frequency grid
--------------

The nominal Fourier spacing is ``1 / T``, where ``T`` is the duration of the
validated time series. ``oversampling`` divides that spacing, while
``nyquist_factor`` controls the upper limit relative to the Nyquist frequency
estimated from the median cadence. The zero-frequency bin is not returned.

Normalization
-------------

``power`` is a one-sided per-bin quantity whose sum equals the variance of the
input flux. ``power_density`` is ``power / frequency_spacing``, so integrating
it over frequency also returns the variance. When flux uncertainties are
provided, inverse-variance weights define both the mean and the variance.

``amplitude`` is ``sqrt(2 * power)`` and therefore has the same units as the
input flux. Oversampled bins are correlated, so individual amplitude values
are not independent Fourier amplitudes when ``oversampling`` is greater than
one.

These definitions deliberately replace the dimensionally inconsistent legacy
amplitude calculation captured by the PBjam compatibility tests.

API
---

.. autofunction:: mimir.power_spectrum

.. autoclass:: mimir.PowerSpectrum
   :members:
