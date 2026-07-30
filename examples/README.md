# Mimir examples

The notebooks are intended to be read in order, but each is self-contained:

1. [`01_local_time_series_and_spectrum.ipynb`](01_local_time_series_and_spectrum.ipynb)
   validates a synthetic light curve, calculates a `nifty-ls` spectrum, and
   checks its normalization.
2. [`02_mast_lightcurve_workflow.ipynb`](02_mast_lightcurve_workflow.ipynb)
   downloads a Kepler light curve and demonstrates both Mimir's target-name
   interface and the individual Lightkurve-backed archive stages.
3. [`03_spectral_window.ipynb`](03_spectral_window.ipynb) compares regular and
   gapped sampling and explains effective independent-frequency spacing.
4. [`04_pbjam_mimir_comparison.ipynb`](04_pbjam_mimir_comparison.ipynb) runs
   one reduced Kepler light curve through both `pbjam.IO.psd` and Mimir,
   comparing the frequency grids, normalization, and power-density spectra.

Install Mimir and the plotting dependency with:

```bash
python -m pip install \
  "mimir-astro @ git+https://github.com/nielsenmb/Mimir.git" matplotlib
```

The MAST and PBjam-comparison notebooks also require optional dependencies:

```bash
python -m pip install \
  "mimir-astro[mast,compat] @ git+https://github.com/nielsenmb/Mimir.git" matplotlib
```

The two synthetic-data notebooks require no network access. The MAST and
PBjam-comparison notebooks download data and may take several minutes depending
on the selected target and products.
