Time-series preparation
=======================

The :class:`mimir.TimeSeries` class is the validated boundary between archive
or user-provided measurements and Mimir's later spectral calculations. It
stores independent NumPy arrays, so callers do not need Lightkurve or a
particular computational backend.

Basic use
---------

.. code-block:: python

   import numpy as np
   from mimir import TimeSeries

   time = np.arange(1000) * 120.0
   flux = np.sin(2 * np.pi * time / 3000.0)

   series = TimeSeries(time, flux, time_unit="s")

   print(series.cadence)
   print(series.duration)
   print(series.duty_cycle)

Cleaning and validation
-----------------------

By default, non-finite samples are removed and samples are sorted by time.
Use ``bad_mask`` to remove known bad measurements, with ``True`` marking a
sample for removal. Set ``invalid="raise"`` or ``sort=False`` when silent
cleaning or sorting is undesirable.

Duplicate timestamps and non-positive finite uncertainties always raise an
error. Combining duplicate measurements requires a scientific choice about
their uncertainties, so Mimir leaves that operation to the caller.

Metadata conventions
--------------------

``cadence`` is the median difference between sorted sample times, and
``duration`` is the difference between the final and initial time. Both use the
unit recorded by ``time_unit``.

``duty_cycle`` compares the retained sample count with the number expected over
the duration at the median cadence. A complete regularly sampled series
therefore has a duty cycle of one. This deliberately corrects the PBjam legacy
endpoint convention, which can return values above one.

API
---

.. autoclass:: mimir.TimeSeries
   :members:
