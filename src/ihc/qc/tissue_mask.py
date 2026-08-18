"""A tissue area that does not depend on which classifier you ran.

Why this module exists
----------------------
Percent area is a ratio, and until 2026-08-17 we were letting the classifier set
*both* halves of it. QuPath omits the ``Ignore*`` class from measurements, so the
tissue area we reported was

    tissue = Abeta area + Negative area

and ``Ignore*`` is a **learned** class. Change the classifier and the denominator
moves with it. Comparing two classifiers on that basis compares two ratios with two
different denominators, which is not a comparison at all.

It was not a small effect on this cohort: 100 of 121 images moved by more than 5%
between the unnormalised and the locally-normalised model, over a range of
0.42-2.40x. Worse, the moving denominator was not symmetric across conditions -- the
candidate classifier gave secondary-only control sections a 21% *larger* denominator
than stained sections from the same animal, in every one of 25 animals. That deflates
control percentages, which is exactly the direction that flatters the conclusion we
were trying to reach. See ADR-0025.

So the denominator is computed here, from DAPI alone, and every classifier is scored
against the same one.

The measure
-----------
Per section, at ~10 um/px (same reasoning as `stain_check` -- below the scale of a
cell, Otsu on DAPI segments nuclei rather than tissue):

1. Otsu on DAPI, with the histogram **clipped at its 99.5th percentile** first.
2. Binary closing (7x7), then fill holes.
3. Keep every connected component larger than 2% of the largest.

Each step earns its place, and two of them were added only after a version without
them failed a validation check. Details in the functions below.

What is established, and what only looks established
----------------------------------------------------
Corrected 2026-08-18 after adversarial review. The original version of this block
presented four numbers as validation. Two of them are consequences of the algorithm
rather than measurements, and one was computed against the wrong channel.

**True by construction, not by measurement.** Otsu after a percentile clip is
*exactly* equivariant under ``x -> a*x + b``: the clip percentile maps affinely, so
the clipped histogram, the optimal threshold and the resulting partition all map
through unchanged. Measured area change under gains of 0.5x-10x and offsets to +150
counts: **0.000%, every time.**

The invariance is in fact broader than affine. Whenever the tissue/glass histogram is
cleanly bimodal, *any* monotonic transform moves the modes and the threshold together
and leaves the partition alone -- verified for sqrt, log and x**1.7, and for
saturation that does not merge the modes. What does move it is anything that erodes
the separation itself: read noise at very low exposure, stray light, saturation that
swallows both modes. Both properties are pinned in ``tests/test_tissue_mask.py``.

Two consequences:

* The repeat-acquisition agreement (0.62% median across a 2.1-5.2x DAPI brightness
  change) is largely *implied*, because those pairs are approximately affine. What
  that test really probes is the non-affine residual -- a real thing, and the mask
  passes, but a much weaker claim than "validated on identical tissue".
* The control/positive ratio of 0.985 is likewise forced. Control sections differ
  from stained ones in the primary antibody, which is a Cy3 effect, and this mask
  reads only DAPI. It could not have come out otherwise.

**Genuinely measured.**

===============================================  =========  ==============
property                                         want       measured
===============================================  =========  ==============
median area vs atlas hemisphere (22.13 mm^2)     close      21.2 mm^2
correlation with **DAPI** brightness, stained    ~0         +0.15
correlation with **DAPI** brightness, controls   ~0         -0.16
tissue area, rapamycin vs control (positives)    ~1.0       0.991 (p=0.26)
===============================================  =========  ==============

The two correlations were previously published as -0.10 and +0.05. Those were
computed against **Cy3** median -- a channel this mask never reads -- because
``fixed_denominator_report.py`` passed the area as the second positional argument
and left ``x`` at its default. Recomputed against DAPI they are larger and of
opposite sign. Neither reaches significance at n = 67 / 52, so the conclusion stands,
but the original check tested the wrong variable.

The last row is new and was never previously checked: a denominator that differed by
treatment arm would shift every percent-area in one direction. It does not, on the
sections that carry the endpoint.

**What remains genuinely load-bearing:** this mask is classifier-independent. Nothing
here reads a classifier output, which is the defect in ``Abeta + Negative`` that
ADR-0025 exists to fix. Claim that, and no more.

Scope
-----
This is whole-section tissue. It is **not** the regional endpoint -- hippocampus and
isocortex come from atlas registration and are a separate mask entirely. This exists
so that classifier comparisons are valid, and so that any percent area we quote has a
denominator whose definition does not change when the numerator does.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .stain_check import CH_DAPI, read_pyramid_plane

__all__ = [
    "TissueArea",
    "OTSU_CLIP_PERCENTILE",
    "CLOSING_SIZE_PX",
    "MIN_COMPONENT_FRACTION",
    "tissue_mask",
    "measure_tissue_area",
]

# Clip the histogram here before thresholding. A single bright speck is enough to
# capture Otsu outright: on one section (tube 32, section 2) the DAPI plane peaks at
# 37,612 counts against a normal 670, and un-clipped Otsu called 0.03% of the frame
# tissue. Otsu maximises between-class variance, so one extreme outlier can make
# "speck vs everything else" the best two-class split available.
OTSU_CLIP_PERCENTILE = 99.5

# Close gaps before taking components. See tissue_mask() -- this is not cosmetic.
CLOSING_SIZE_PX = 7

# Keep every component above this fraction of the largest. Not "keep the largest":
# sections legitimately arrive in two or more pieces.
MIN_COMPONENT_FRACTION = 0.02


@dataclass
class TissueArea:
    """One section's tissue area, and how it was obtained."""

    ets_path: str
    level: int = -1
    um_per_px: float = float("nan")
    area_mm2: float = float("nan")
    n_components: int = 0
    problems: list[str] = field(default_factory=list)


def tissue_mask(plane: np.ndarray) -> np.ndarray:
    """Boolean tissue mask from one DAPI plane.

    ``plane`` carries NaN at never-acquired tile positions (see
    `read_pyramid_plane`); those are never tissue, and NaN-aware statistics are used
    throughout so they cannot drag the threshold.

    Two of these steps exist because a simpler version failed validation, and both
    failures are worth stating because both looked fine until measured.

    **Closing and multi-component selection.** Otsu alone, keeping only the largest
    connected component, gave a mask whose area correlated with section brightness at
    +0.38 (stained) and +0.43 (controls). The threshold was not the problem -- the
    mask was. Dim sections fragment into many small pieces at any threshold, so
    "largest component" discards more area the dimmer the section, and the brightness
    dependence was manufactured by the cleanup step rather than measured from the
    tissue. Closing bridges the fragments; the 2% rule then keeps genuine separate
    pieces without readmitting debris.

    **Otsu rather than a per-section z-score.** The obvious scale-invariant choice --
    threshold at ``median + 3*MAD`` of the section's own DAPI -- collapses to a median
    of 0.6 mm^2 against a true ~21 mm^2. Tissue occupies roughly 55% of these frames,
    so the global median sits *inside* tissue and a z-threshold selects only the
    brightest nuclei. Invariance by construction buys nothing if the statistic is
    computed over the wrong population.
    """
    from scipy import ndimage as ndi
    from skimage.filters import threshold_otsu

    finite = plane[np.isfinite(plane)]
    if finite.size == 0:
        return np.zeros(plane.shape, dtype=bool)

    ceiling = np.percentile(finite, OTSU_CLIP_PERCENTILE)
    threshold = threshold_otsu(np.minimum(finite, ceiling))

    mask = np.isfinite(plane) & (plane > threshold)
    mask = ndi.binary_closing(mask, structure=np.ones((CLOSING_SIZE_PX, CLOSING_SIZE_PX)))
    mask = ndi.binary_fill_holes(mask)

    labels, n = ndi.label(mask)
    if n == 0:
        return mask
    sizes = ndi.sum(mask, labels, range(1, n + 1))
    keep = np.where(sizes > MIN_COMPONENT_FRACTION * sizes.max())[0] + 1
    mask = np.isin(labels, keep)

    # Re-intersect with acquired support. `binary_fill_holes` will happily fill an
    # interior block of never-acquired tiles and hand it back as tissue -- measured at
    # 400 of 400 pixels on a synthetic 20x20 NaN block before this line existed. The
    # scanner's sample mask skips tile positions with no tissue, so those pixels are
    # missing support, not dark tissue, and they must not enter the denominator
    # (ADR-0010, ADR-0024).
    return mask & np.isfinite(plane)


def measure_tissue_area(
    ets_path: str,
    *,
    pixel_size_um: float,
    true_size: tuple[int, int] | None = None,
) -> TissueArea:
    """Tissue area in mm^2 for one section, from DAPI only.

    Args:
        pixel_size_um: base-level calibration, from the manifest.
        true_size: base-level ``(width_px, height_px)`` from tag 2053. Pass it -- the
            tile grid is padded with zeros that are indistinguishable from black
            tissue, and the pad is large at coarse levels.
    """
    result = TissueArea(ets_path=ets_path)
    try:
        plane, level = read_pyramid_plane(
            ets_path, CH_DAPI, pixel_size_um=pixel_size_um, true_size=true_size
        )
    except ValueError as exc:
        result.problems.append(str(exc))
        return result

    result.level = level
    result.um_per_px = float(pixel_size_um) * (2 ** level)

    mask = tissue_mask(plane)
    if not mask.any():
        result.problems.append("no tissue found")
        return result

    from scipy import ndimage as ndi

    _, result.n_components = ndi.label(mask)
    result.area_mm2 = float(mask.sum()) * (result.um_per_px ** 2) / 1e6

    # A section that comes back far from the cohort's ~21 mm^2 is a read to inspect,
    # not a number to trust. Flag rather than raise: the caller aggregates.
    if result.area_mm2 < 5.0:
        result.problems.append(f"implausibly small tissue area ({result.area_mm2:.2f} mm2)")

    return result
