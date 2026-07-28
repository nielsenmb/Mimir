"""Tests for spectral-window and effective-spacing calculations."""

import numpy as np
import pytest

from mimir.window import (
    SpectralWindow,
    effective_frequency_spacing,
    spectral_window,
)


def test_regular_sampling_window_is_normalized_and_symmetric():
    """A regular sampling pattern should produce a symmetric unit-height window."""
    time = np.arange(512) * 120.0 / 86400.0

    result = spectral_window(time, oversampling=8)

    assert isinstance(result, SpectralWindow)
    assert result.power[result.n_bins // 2] == pytest.approx(1.0)
    assert result.frequency[result.n_bins // 2] == pytest.approx(0.0)
    assert result.frequency == pytest.approx(-result.frequency[::-1])
    assert result.power == pytest.approx(result.power[::-1], abs=1e-14)
    assert result.frequency_unit == "uHz"


def test_regular_sampling_effective_spacing_approaches_nominal_spacing():
    """A long uninterrupted series should recover its nominal Fourier spacing."""
    time = np.arange(1024) * 120.0 / 86400.0

    result = spectral_window(time, oversampling=10)

    assert result.effective_frequency_spacing == pytest.approx(
        result.nominal_frequency_spacing,
        rel=5e-3,
    )
    assert effective_frequency_spacing(time) == pytest.approx(
        result.effective_frequency_spacing
    )


def test_gaps_increase_effective_spacing():
    """Removing observations should broaden the integrated sampling window."""
    complete = np.arange(1024) * 120.0 / 86400.0
    sample = np.arange(complete.size)
    gapped = complete[~((sample > 400) & (sample < 700))]

    complete_spacing = effective_frequency_spacing(complete)
    gapped_spacing = effective_frequency_spacing(gapped)

    assert gapped_spacing > complete_spacing


def test_window_is_invariant_to_time_translation():
    """Shifting every timestamp should not change the spectral-window power."""
    time = np.arange(256) * 120.0 / 86400.0

    original = spectral_window(time, oversampling=4)
    shifted = spectral_window(time + 2_450_000.0, oversampling=4)

    assert shifted.frequency == pytest.approx(original.frequency)
    assert shifted.power == pytest.approx(original.power, abs=1e-9)
    assert shifted.effective_frequency_spacing == pytest.approx(
        original.effective_frequency_spacing,
        rel=1e-8,
    )


@pytest.mark.parametrize("oversampling", [True, 0, 1.5])
def test_invalid_oversampling_is_rejected(oversampling):
    """Oversampling must be a positive integer."""
    time = np.arange(32, dtype=float)

    expected = TypeError if oversampling is True or oversampling == 1.5 else ValueError
    with pytest.raises(expected):
        spectral_window(time, oversampling=oversampling)


def test_too_narrow_window_is_rejected():
    """The requested half-width must include bins either side of zero."""
    time = np.arange(32, dtype=float)

    with pytest.raises(ValueError, match="at least one frequency bin"):
        spectral_window(time, half_width=1e-9)
