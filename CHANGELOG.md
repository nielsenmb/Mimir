# Changelog

All notable changes to Mimir are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/) and the
project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Concise target calls such as `power_spectrum(target, mast_kwargs)` with flat
  Lightkurve search filters and separate optional `lightcurve_kwargs`.
- Flat filter mappings for `search_lightcurves`, with documented product
  inspection and selection before downloading.
- Direct MAST target-name inputs for time-series, power-spectrum, and
  spectral-window workflows.
- Explicit, regularly spaced frequency grids for power-spectrum calculation,
  allowing multiple targets to be evaluated on identical frequency bins.

### Changed

- Retained positional arrays, distinct named input forms, and legacy nested
  MAST options while simplifying the target workflow.
- Split numerical inputs into explicit, mutually exclusive ``time``/``flux``,
  ``time_series``, and ``target`` arguments.

### Fixed

- Fixed PSD density normalization to use a reference grid independent of
  output spacing, oversampling, and frequency range. Only the default automatic
  reference grid is required to integrate exactly to the chosen variance.
- Preserve NumPy and Astropy masks in array, sampling-window, and Lightkurve
  inputs; exclude masked samples before estimating the ppm normalization.
- Convert Lightkurve time coordinates and cadence using Astropy time arithmetic
  and physical units, independent of JD, mission, Unix, or calendar display
  formats. Absolute times now use elapsed time since JD 2451545.0 in the input
  scale; unitless numeric columns remain interpreted as days.
- Named Lightkurve exposure selectors (`short`, `long`, `fast`) no longer enter
  numerical exposure-time validation during default reduction.

## [0.1.0] - 2026-07-28

### Added

- Validated NumPy-based `TimeSeries` containers with masking, cleaning,
  sorting, uncertainty validation, and cadence metadata.
- One-sided Lomb–Scargle power spectra calculated directly with nifty-ls.
- Parseval-normalized power, power density, and sinusoidal semi-amplitude.
- Spectral-window and effective independent-frequency-spacing calculations.
- Optional Lightkurve-backed MAST search, download, reduction, and conversion.
- PBjam 2.0.4 characterisation tests recording the extraction baseline.
- Test, lint, build, and Sphinx documentation workflows for Python 3.10–3.13.

### Changed

- Corrected the legacy PBjam duty-cycle endpoint convention.
- Corrected weighted centring and the dimensionally inconsistent legacy
  amplitude definition.
- Made spectrum normalization independent of oversampling and requested
  super-Nyquist extension.

[Unreleased]: https://github.com/nielsenmb/Mimir/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/nielsenmb/Mimir/releases/tag/v0.1.0
