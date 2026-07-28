PBjam compatibility contract
============================

Mimir was extracted from ``pbjam.IO`` by first recording the behaviour of
PBjam 2.0.4 and then implementing the corresponding features behind a new,
validated API. The tests in ``tests/compatibility`` remain characterisation
tests: they document the migration baseline rather than prescribing Mimir's
final numerical behaviour.

Recorded legacy behaviour
-------------------------

The compatibility suite records:

* removal of non-finite and user-masked samples;
* cadence, duration, sample count, and legacy duty-cycle calculation;
* the Lomb--Scargle frequency grid and Nyquist truncation;
* unweighted Parseval normalization for power and power density;
* recovery of the dominant frequency of a synthetic sinusoid;
* the legacy weighted normalization and amplitude conventions;
* time-domain window construction and its duplicated padding endpoint; and
* mocked Lightkurve search, download, and stitching calls.

Intentional Mimir differences
-----------------------------

Mimir does not reproduce known defects merely for compatibility. Its current
API:

* bounds duty cycles at one;
* uses inverse-variance-weighted centring when uncertainties are present;
* reports a dimensionally consistent sinusoidal semi-amplitude;
* keeps amplitude stable under frequency-grid oversampling;
* normalizes using the physical one-sided band through Nyquist; and
* computes the sampling window directly from the Fourier transform of the
  observation mask.

These differences are covered by Mimir's ordinary regression tests. The
Astropy high-frequency bias is also not a compatibility target: Mimir calls
``nifty-ls`` directly.

Running the tests
-----------------

The PBjam dependency is isolated in the ``compat`` development extra:

.. code-block:: console

   python -m pip install ".[test,compat]"
   pytest -m compatibility

Ordinary Mimir installations and tests do not require PBjam. Continuous
integration runs the compatibility suite separately on Python 3.12.
