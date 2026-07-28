PBjam compatibility contract
============================

Mimir is being extracted from ``pbjam.IO`` in two distinct phases:

#. record the behaviour of the existing implementation;
#. implement the equivalent Mimir feature and compare it with that record.

The tests in ``tests/compatibility`` are therefore characterisation tests.
They describe PBjam 2.0.4 rather than prescribing Mimir's final behaviour.
This distinction matters because the inherited implementation contains known
numerical defects. Reproducing a result during extraction and endorsing that
result scientifically are not the same thing.

Covered behaviour
-----------------

The initial contract records:

* removal of non-finite and user-masked samples;
* cadence, duration, sample count, and legacy duty-cycle calculation;
* the Lomb--Scargle frequency grid and Nyquist truncation;
* unweighted Parseval normalization for power and power density;
* recovery of the dominant frequency of a synthetic sinusoid;
* the present weighted normalization convention;
* the present amplitude calculation;
* time-domain window construction and its legacy duplicated padding endpoint;
  and
* mocked Lightkurve search, download, and stitching calls.

The amplitude calculation, the weighted centring convention, and duty cycles
above one are intentionally captured as legacy behaviour. They will be
corrected only in the later numerical-audit stage, with explicit regression
tests and release notes.

Running the tests
-----------------

The PBjam dependency is isolated in the ``compat`` development extra:

.. code-block:: console

   python -m pip install ".[test,compat]"
   pytest -m compatibility

Ordinary Mimir tests do not require PBjam. Continuous integration runs the
compatibility suite separately on Python 3.12 so the migration boundary
remains visible.

The Astropy high-frequency bias is not a compatibility target. Mimir will use
``nifty-ls`` when the spectrum implementation is introduced.
