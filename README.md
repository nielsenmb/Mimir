# Mimir

Mimir provides shared tools for preparing astronomical time series and
computing power-density spectra for asteroseismology. It is being extracted
from `pbjam.IO` so that PBjam, Skuld, Urdr, and other projects can use the same
well-tested implementation.

The package is currently under development. It provides a validated
`TimeSeries` container and a Parseval-normalized power-spectrum function:

```python
from mimir import TimeSeries, power_spectrum

series = TimeSeries(time, flux, flux_err, time_unit="d")
spectrum = power_spectrum(series, oversampling=2)
print(spectrum.frequency, spectrum.power_density)
```

MAST access and basic Lightkurve reduction are available through the optional
`mast` dependency:

```python
from mimir import load_lightcurve

series = load_lightcurve(
    "KIC 8006161",
    search_kwargs={"mission": "Kepler", "exptime": 60},
    numax=3500,
)
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

Time-series validation, nifty-ls power-spectrum calculation, and optional
Lightkurve-based MAST access and reduction are implemented.

## License

Mimir is distributed under the MIT License.
