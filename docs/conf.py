"""Sphinx configuration for Mimir."""

from importlib.metadata import version as get_version

project = "Mimir"
author = "Martin Nielsen"
copyright = "2026, Martin Nielsen"
version = release = get_version("mimir-astro")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

html_theme = "sphinx_rtd_theme"
exclude_patterns = ["_build"]
