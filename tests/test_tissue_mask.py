"""Tests for the classifier-independent tissue denominator.

Three of these pin defects that a working-looking earlier version actually had, and
one pins the property the whole module exists for. In order of what they cost us:

- `test_area_is_invariant_to_section_brightness` -- the point. A denominator that
  scales with brightness reintroduces exactly the artefact it was built to remove.
- `test_fragmented_dim_section_keeps_its_area` -- the +0.38 correlation in ADR-0025
  came from cleanup, not from thresholding: keeping only the largest component
  discards more area the dimmer the section.
- `test_bright_speck_does_not_capture_otsu` -- one 37,612-count dust speck against a
  normal 670 was enough to make Otsu call 0.03% of a frame tissue.
- `test_two_genuine_pieces_are_both_kept` -- sections do arrive in pieces.
"""

from __future__ import annotations

import numpy as np
import pytest

from ihc.qc.tissue_mask import tissue_mask


def _section(brightness: float = 200.0, size: int = 200, glass: float = 5.0) -> np.ndarray:
    """Tissue occupying the middle ~55% of the frame, on glass."""
    plane = np.full((size, size), glass, dtype=np.float32)
    lo, hi = int(size * 0.22), int(size * 0.78)
    plane[lo:hi, lo:hi] = brightness
    return plane


def test_area_is_invariant_to_section_brightness():
    """A 10x gain change must not change the measured tissue area."""
    dim = tissue_mask(_section(brightness=60.0, glass=1.5))
    bright = tissue_mask(_section(brightness=600.0, glass=15.0))
    assert dim.sum() == bright.sum()


def test_fragmented_dim_section_keeps_its_area():
    """Closing must bridge gaps that a dim section breaks into.

    Without it the pipeline took the largest component and silently shed area, which
    is how a mask-fragility artefact was mistaken for a brightness effect.
    """
    plane = _section(brightness=80.0)
    plane[::6, :] = 5.0  # thin dropouts, as a dim section fragments
    mask = tissue_mask(plane)
    whole = tissue_mask(_section(brightness=80.0))
    assert mask.sum() > 0.9 * whole.sum()


def test_bright_speck_does_not_capture_otsu():
    """One extreme outlier must not become the best two-class split."""
    plane = _section(brightness=200.0)
    clean = tissue_mask(plane)
    plane[3, 3] = 37_612.0  # the real value from tube 32, section 2
    assert tissue_mask(plane).sum() == pytest.approx(clean.sum(), rel=0.02)


def test_two_genuine_pieces_are_both_kept():
    plane = np.full((200, 200), 5.0, dtype=np.float32)
    plane[20:90, 20:180] = 200.0
    plane[110:180, 20:180] = 200.0
    mask = tissue_mask(plane)
    from scipy import ndimage as ndi

    assert ndi.label(mask)[1] == 2
    assert mask.sum() > 0.9 * (70 * 160 * 2)


def test_never_acquired_support_is_not_tissue():
    """NaN marks tile positions the scanner skipped; they are never tissue."""
    plane = _section(brightness=200.0)
    plane[:40] = np.nan
    mask = tissue_mask(plane)
    assert not mask[:40].any()


def test_empty_plane_returns_empty_mask():
    plane = np.full((50, 50), np.nan, dtype=np.float32)
    assert tissue_mask(plane).sum() == 0


def test_mask_is_exactly_invariant_to_gain_and_offset():
    """Pins the property that makes the repeat-pair agreement *not* independent evidence.

    Otsu after a percentile clip is exactly equivariant under ``x -> a*x + b``: the clip
    percentile maps affinely, so the clipped histogram, the threshold and the partition
    all map through unchanged. This is a genuine strength -- the denominator cannot drift
    with exposure -- but it also means a repeat-acquisition test on affinely-related
    images cannot fail, and must not be quoted as validation. Documented in the module
    docstring; pinned here so nobody re-derives that argument from a lucky measurement.
    """
    base = _section(brightness=200.0)
    n = tissue_mask(base).sum()
    for a, b in [(0.5, 0.0), (2.0, 0.0), (10.0, 0.0), (1.0, 100.0), (6.0, 120.0)]:
        assert tissue_mask((base * a + b).astype(np.float32)).sum() == n


def test_invariance_is_broader_than_affine_on_a_bimodal_histogram():
    """Any monotonic transform, not just affine, leaves a clean two-mode mask alone.

    Otsu picks a threshold between the modes; a monotonic map moves the modes and the
    threshold together, so the partition survives. This is why the repeat-acquisition
    agreement has even less evidential value than the affine argument alone implies.
    """
    base = _section(brightness=200.0)
    n = tissue_mask(base).sum()
    for f in (np.sqrt, np.log1p, lambda x: x ** 1.7):
        assert tissue_mask((f(np.maximum(base, 0)) * 20.0).astype(np.float32)).sum() == n


def test_what_actually_breaks_the_mask_is_noise_not_scaling():
    """Rescaling cannot move this mask; degrading the separation can.

    Worth pinning because it bounds what a repeat-acquisition test can ever detect. Two
    acquisitions that differ in exposure produce the same mask by construction. Only a
    change that erodes the tissue/glass separation -- read noise at very low exposure,
    stray light, saturation merging the modes -- moves it. So the mask's stability across
    the rescan pairs is a statement about the algorithm, not about the acquisitions.
    """
    rng = np.random.default_rng(0)
    base = _section(brightness=200.0)          # glass 5, tissue 200: cleanly separated
    n = tissue_mask(base).sum()
    mild = base + rng.normal(0, 5, base.shape)
    assert tissue_mask(mild.astype(np.float32)).sum() == pytest.approx(n, rel=0.02)
    swamped = base + rng.normal(0, 400, base.shape)   # modes no longer separable
    assert tissue_mask(swamped.astype(np.float32)).sum() != pytest.approx(n, rel=0.05)
