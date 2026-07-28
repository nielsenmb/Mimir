# Mimir

Mimir provides shared tools for preparing astronomical time series and
computing power-density spectra for asteroseismology. It is being extracted
from `pbjam.IO` so that PBjam, Skuld, Urdr, and other projects can use the same
well-tested implementation.

The package is currently under development. Its first scientific component is
the validated `TimeSeries` container:

```python
from mimir import TimeSeries

series = TimeSeries(time, flux, flux_err, time_unit="d")
print(series.cadence, series.duration, series.duty_cycle)
```

## Design

Mimir separates three responsibilities:

- validation and preparation of time-series data;
- accurate Lomb–Scargle power-spectrum calculation using
  [nifty-ls](https://github.com/flatironinstitute/nifty-ls);
- optional MAST access and light-curve reduction through Lightkurve.

NumPy arrays are the public interchange format. Lightkurve is optional so that
local spectrum calculations do not require archive-access dependencies.

## Development installation

Clone the repository and install the development dependencies:

```bash
python -m pip install -e ".[test,mast,docs]"
```

Run the checks with:

```bash
ruff check .
pytest
python -m build
```

The distribution name is `mimir-astro`, because the unrelated name `mimir` is
already registered on PyPI. The Python import remains:

```python
import mimir
```

## Status

Time-series validation and metadata are implemented. Power-spectrum
calculation and archive access will be introduced in subsequent, independently
tested changes.

## License

Mimir is distributed under the MIT License.
