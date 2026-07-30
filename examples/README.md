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

Install Mimir and the plotting dependency with:

```bash
python -m pip install \
  "mimir-astro @ git+https://github.com/nielsenmb/Mimir.git" matplotlib
```

The MAST notebook also requires the optional archive dependency:

```bash
python -m pip install \
  "mimir-astro[mast] @ git+https://github.com/nielsenmb/Mimir.git" matplotlib
```

The two synthetic-data notebooks require no network access. The MAST notebook
downloads data and may take several minutes depending on the selected target
and products.
