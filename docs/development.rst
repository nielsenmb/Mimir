Development plan
================

The implementation will be developed in small, testable stages:

#. Characterise the behaviour of ``pbjam.IO``.
#. Add validated time-series containers.
#. Add nifty-ls power-spectrum calculations and normalization.
#. Add optional Lightkurve-based MAST access and reduction.
#. Audit numerical behaviour before changing PBjam to depend on Mimir.

Compatibility tests will establish that moving the code does not accidentally
change scientific results. Deliberate bug fixes and numerical changes will
then be introduced with separate regression tests.
