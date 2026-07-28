"""Tests for the initial package structure."""

from importlib import import_module
from importlib.metadata import version

import mimir


def test_version_matches_distribution_metadata():
    """The public version should match the installed distribution."""
    assert mimir.__version__ == version("mimir-astro")


def test_planned_modules_are_importable():
    """The initial scientific module boundaries should remain importable."""
    for module in ("mast", "spectrum", "timeseries"):
        assert import_module(f"mimir.{module}") is not None
