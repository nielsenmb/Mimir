# Changelog

All notable changes to Mimir are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/) and the
project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Direct MAST target-name inputs for time-series, power-spectrum, and
  spectral-window workflows.

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
