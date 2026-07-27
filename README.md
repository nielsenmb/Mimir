# Mimir

Mimir will provide shared tools for preparing astronomical time series and
computing power-density spectra for asteroseismology. It is being extracted
from `pbjam.IO` so that PBjam, Skuld, Urdr, and other projects can use the same
well-tested implementation.

The package is currently under development and does not yet provide a stable
scientific API.

## Design

Mimir will separate three responsibilities:

- validation and preparation of time-series data;
- accurate Lomb–Scargle power-spectrum calculation using
  [nifty-ls](https://github.com/flatironinstitute/nifty-ls);
- optional MAST access and light-curve reduction through Lightkurve.

NumPy arrays will be the public interchange format. Lightkurve is optional so
that local spectrum calculations do not require archive-access dependencies.

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

The eventual distribution name is `mimir-astro`, because the unrelated name
`mimir` is already registered on PyPI. The Python import remains:

```python
import mimir
```

## Status

The first development stage establishes the package, documentation, tests, and
continuous integration. Scientific functionality will be introduced in
subsequent, independently tested changes.

## License

Mimir is distributed under the MIT License.
