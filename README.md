# Mimir

[![Tests](https://github.com/nielsenmb/Mimir/actions/workflows/tests.yml/badge.svg)](https://github.com/nielsenmb/Mimir/actions/workflows/tests.yml)
[![Python 3.10–3.13](https://img.shields.io/badge/python-3.10--3.13-blue.svg)](https://www.python.org/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Mimir provides shared tools for preparing astronomical time series and
computing power-density spectra for asteroseismology. It separates numerical
spectrum calculation from optional archive access so that projects such as
PBjam, Skuld, and Urdr can use the same tested implementation.

Mimir currently provides:

- validated NumPy-based time-series containers;
- Lomb–Scargle power spectra calculated directly with
  [nifty-ls](https://github.com/flatironinstitute/nifty-ls);
- explicit Parseval normalization, power density, and amplitude;
- spectral-window and effective independent-frequency-spacing calculations;
- optional MAST access and basic reduction through Lightkurve.

## Installation

Clone the repo and do:

```bash
python -m pip install -e .
```
for the light-weight version, or for if you want to include download handling from MAST :

```bash
python -m pip install -e ".[mast]"
```
The import the module:

```python
import mimir
```

Until the first PyPI release, install the latest public source with:

```bash
python -m pip install "mimir-astro @ git+https://github.com/nielsenmb/Mimir.git"
```

## Quick start

```python
from mimir import TimeSeries, power_spectrum, spectral_window

series = TimeSeries(time, flux, flux_err, time_unit="d", flux_unit="ppm")
spectrum = power_spectrum(time_series=series, oversampling=2)
window = spectral_window(time_series=series)

print(spectrum.frequency, spectrum.power_density)
print(window.effective_frequency_spacing)
```

Use an explicit regular grid when spectra for several targets must be directly
comparable bin by bin:

```python
shared_frequency = np.arange(10.0, 5000.0, 0.1)  # uHz
spectrum = power_spectrum(time_series=series, frequency=shared_frequency)
```

MAST access and basic Lightkurve reduction are available with the `mast`
extra:

```python
from mimir import power_spectrum, search_lightcurves

target = "KIC 8006161"
mast_kwargs = {"mission": "Kepler", "author": "Kepler", "exptime": 60}
spectrum = power_spectrum(target, mast_kwargs)
```

To inspect available light-curve products before downloading anything:

```python
products = search_lightcurves(target)
print(products)  # Available authors/pipelines, observing segments, and exposure times
products.table  # Full metadata in a notebook
```

`search_lightcurves` searches the products supported by Lightkurve, including
supported high-level science light curves. Use its returned metadata to choose
filters such as `author`, `mission`, `exptime`, `quarter`, or `sector`. The
quick PSD call downloads all matching products; choose a consistent pipeline
and cadence when different products cover the same observations. The
[MAST guide](docs/mast.rst) also shows how to select individual result rows.

`mast_kwargs` is a **flat dictionary of search filters**. Optional reduction
settings are separate:

```python
spectrum = power_spectrum(target, mast_kwargs, lightcurve_kwargs={"numax": 3500})
```

Named `target=`, `time=`/`flux=`, and `time_series=` inputs remain available and
mutually exclusive. Positional arrays remain supported. Earlier nested
`mast_kwargs={"search_kwargs": {...}, "numax": ...}` calls still work, but are
not needed for new code.

## Example notebooks

The [`examples`](examples) directory contains tutorials for:

- validating local arrays and calculating a normalized power spectrum;
- downloading and reducing a MAST light curve through Lightkurve; and
- interpreting spectral windows and effective frequency spacing for gapped
  observations; and
- comparing the PBjam and Mimir spectra of the same reduced Kepler target.

The local-array and spectral-window notebooks use deterministic synthetic data
and run without network access.

## Numerical conventions

Power is one-sided and normalized over the physical band through Nyquist. At
the default `nyquist_factor=1`, the sum of the power equals the input flux
variance. Power density integrates to the same variance, and amplitude is a
sinusoidal semi-amplitude with the same units as the input flux. See
[the spectrum documentation](docs/spectrum.rst) for the complete definitions.

Mimir uses nifty-ls directly rather than Astropy's standard Lomb–Scargle
backend. NumPy arrays are the public interchange format; Lightkurve is required
only for archive access and reduction.

## Contributing

Bug reports, numerical test cases, and focused pull requests are welcome. See
[CONTRIBUTING.md](CONTRIBUTING.md) for the development setup and conventions.
All public functions and classes use NumPy-style docstrings.

## Citation

If Mimir contributes to published research, please cite the software using the
metadata in [CITATION.cff](CITATION.cff). Release-specific archive information
can be added after the first public release.

## License

Mimir is distributed under the [MIT License](LICENSE).
