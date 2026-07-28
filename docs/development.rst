Development
===========

Architecture
------------

Mimir keeps its main responsibilities separate:

* :mod:`mimir.timeseries` validates and prepares NumPy-based measurements;
* :mod:`mimir.spectrum` computes normalized spectra with nifty-ls;
* :mod:`mimir.window` describes the sampling pattern and effective frequency
  spacing; and
* :mod:`mimir.mast` provides optional Lightkurve-backed archive access and
  reduction.

This boundary keeps local numerical work independent of Lightkurve and lets
downstream packages choose their own retrieval and reduction defaults.

Development checks
------------------

Install the complete development environment and run:

.. code-block:: console

   python -m pip install -e ".[test,mast,docs]"
   ruff check .
   pytest
   python -m build
   sphinx-build -W -b html docs docs/_build/html

The PBjam characterisation suite is isolated in a separate optional
environment:

.. code-block:: console

   python -m pip install -e ".[test,compat]"
   pytest -m compatibility

Contribution conventions
------------------------

All public functions and classes use NumPy-style docstrings. Numerical changes
must include regression tests and document any effect on units, normalization,
or downstream thresholds. Ordinary tests must not make live MAST requests.

See ``CONTRIBUTING.md`` for the complete contribution guidance.
