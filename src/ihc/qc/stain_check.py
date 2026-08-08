"""Cross-check the recorded staining condition against the pixels.

`condition` -- positive or negative -- is the one manifest column whose corruption
produces silently wrong science, and it comes from a hand-transcribed bench record
(`config/slides.csv`). Nothing downstream can detect a transcription error in it: a
positive section mislabelled negative simply becomes a very bright control.

So this module measures, per section, how much Aß-channel signal is actually there,
and reports whether that agrees with the record.

**It does not, and must not, assign `condition`.** Deriving the staining condition
from the pixels is the forbidden inference: a section whose primary stain genuinely
failed would be silently relabelled a control, converting a bench failure into a
clean-looking data point. This is a cross-check that raises a disagreement for a
human to settle -- nothing more. See ADR-0009.

Use the GFAP index, not the Aß index
------------------------------------
Both are computed, but **`gfap_index` is the discriminator**: measured over this
cohort it separates recorded positives from recorded negatives perfectly
(AUC 1.000, 0 of 128 misassigned, medians 9.6x apart), while `abeta_index` does not
(AUC 0.966, 7 of 128 misassigned).

That is not a quirk of the arithmetic, it is the biology. The negative control here
is secondary-antibody-only, so *both* markers drop together -- which means either
channel could in principle serve. But GFAP is abundant and spread through every
section, whereas Aß is sparse, focal, and genuinely varies between animals by
disease burden and by age. A low Aß index is therefore ambiguous: it may mean the
primary antibody was absent, or it may mean this animal simply has few plaques.
A low GFAP index has no such competing explanation.

Read `abeta_index` as context, not as evidence.

The measure
-----------
Per section, at a coarse pyramid level (this is a gross signal check, not
morphometry):

    index = (marker signal per ms) / (DAPI signal per ms)

Cy3 signal is the 99.5th percentile inside tissue minus the glass background --
a high percentile because plaques are sparse and bright, so a median would be
dominated by unstained parenchyma. DAPI signal is the *median* inside tissue,
which tracks how much tissue and how well it took stain at all.

Both are divided by their exposure time because exposure varies up to 12.6x across
this cohort (ADR-0004). Dividing by DAPI then removes what is common to both
channels -- section thickness, focus, lamp drift -- leaving something comparable
between animals.

The separation between positives and negatives is an *observed property of this
cohort*, not a designed threshold. Report it; do not hard-code it. An earlier
ad-hoc version of this check measured on the padded tile grid and reported a
separation figure that could not be reproduced once the padding was cropped
(ADR-0020) -- which is why this now lives in a module with tests instead of in a
terminal.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from ..ingest.verify import read_ets_index

__all__ = [
    "SectionStain",
    "MIN_PLANE_WIDTH_PX",
    "read_pyramid_plane",
    "measure_section",
    "audit_manifest",
]

# Channel order in this cohort, verified visually (ADR-0003).
CH_DAPI, CH_FITC, CH_CY3 = 0, 1, 2

# Measure at a fixed *physical* scale, not a fixed pixel size, so the number means
# the same thing on a large section and a small one.
#
# ~10 um/px is chosen deliberately and the choice is load-bearing. The tissue mask
# comes from Otsu on DAPI, and DAPI at native resolution is discrete bright nuclei on
# dark neuropil -- so at 0.325 um/px Otsu segments *nuclei*, not tissue, and the
# apparent tissue fraction collapses (63% at 10 um/px against 4% at 1.3 um/px on the
# same section). Only once blurred to a scale coarser than a cell does "above
# threshold" mean "tissue". This is a gross signal check, not morphometry, so
# measuring coarse is right -- but it is right for a reason, not for speed.
TARGET_UM_PER_PX = 10.0

# Percentiles: high for the sparse bright marker, median for the counterstain.
CY3_PERCENTILE = 99.5
BACKGROUND_PERCENTILE = 5.0


@dataclass
class SectionStain:
    """What one section's pixels say about its staining."""

    ets_path: str
    level: int = -1
    plane_width_px: int = 0
    abeta_index: float = float("nan")
    gfap_index: float = float("nan")
    tissue_fraction: float = float("nan")
    problems: list[str] = field(default_factory=list)


def _decode_tile(raw: bytes) -> np.ndarray | None:
    """Decode one JPEG2000 tile, or return None if it cannot be decoded.

    A tile that fails to decode is a data point about the file, not a reason to
    abandon the section -- so the caller counts them rather than raising.
    """
    try:
        import imagecodecs
    except ImportError as exc:  # pragma: no cover - doctor checks this
        raise RuntimeError("imagecodecs is required to read pixels") from exc
    try:
        arr = imagecodecs.jpeg2k_decode(raw)
    except Exception:
        return None
    if arr is None:
        return None
    if arr.ndim == 3:
        arr = arr[..., 0]
    return arr


def read_pyramid_plane(
    ets_path: str,
    channel: int,
    *,
    level: int | None = None,
    pixel_size_um: float | None = None,
    true_size: tuple[int, int] | None = None,
) -> tuple[np.ndarray, int]:
    """Assemble one channel of one pyramid level into a single array.

    Tile positions that were never acquired -- 1.4-23.6% of the bounding box in this
    cohort -- come back as NaN, **not** zero. That distinction is the whole point:
    zero would be read as legitimately black tissue and would drag every percentile
    down. Callers must use nan-aware statistics.

    Args:
        true_size: ``(width_px, height_px)`` of the **base** level, from tag 2053 via
            the manifest. Pass it. The tile grid is padded out to a whole number of
            512 px tiles, and the pad is stored as zeros indistinguishable from real
            black pixels -- at level 6 a 17,920 px image occupies 280 px of a single
            512 px tile, so more than two thirds of the assembled plane is pad. Left
            uncropped, that pad sits below every real pixel and drags the background
            percentile to zero, which silently inflates every signal measured against
            it. This is the same padded-dimensions defect `verify.py` documents.

    Returns:
        ``(plane, level)`` -- the assembled float32 array and the level used.

    Raises:
        ValueError: the requested channel or level holds no tiles.
    """
    index = read_ets_index(ets_path)
    chunks = index["chunks"]
    if not chunks:
        raise ValueError(f"no tiles in {os.path.basename(ets_path)}")

    n_dim = index["n_dim"]
    ch_axis = 2 if n_dim >= 4 else None
    tile_w, tile_h = index["tile"][0], index["tile"][1]
    levels = sorted({c[0][-1] for c in chunks})

    def tiles_at(lv: int) -> list[tuple[tuple[int, ...], int, int]]:
        return [
            c for c in chunks
            if c[0][-1] == lv and (ch_axis is None or c[0][ch_axis] == channel)
        ]

    if level is None:
        if not pixel_size_um:
            raise ValueError("pass either level or pixel_size_um")
        # Level whose scale is nearest TARGET_UM_PER_PX, clamped to what exists.
        want = max(0, round(np.log2(TARGET_UM_PER_PX / float(pixel_size_um))))
        available = [lv for lv in levels if tiles_at(lv)]
        if not available:
            raise ValueError(f"channel {channel} has no tiles in {ets_path}")
        level = min(available, key=lambda lv: (abs(lv - want), lv))

    sel = tiles_at(level)
    if not sel:
        raise ValueError(f"channel {channel} has no tiles at level {level} in {ets_path}")

    cols = max(c[0][0] for c in sel) + 1
    rows = max(c[0][1] for c in sel) + 1
    plane = np.full((rows * tile_h, cols * tile_w), np.nan, dtype=np.float32)

    with open(ets_path, "rb") as fh:
        for coords, offset, nbytes in sel:
            if offset < 0 or nbytes <= 0 or offset + nbytes > index["file_bytes"]:
                continue
            fh.seek(offset)
            tile = _decode_tile(fh.read(nbytes))
            if tile is None:
                continue
            y0, x0 = coords[1] * tile_h, coords[0] * tile_w
            h = min(tile.shape[0], plane.shape[0] - y0)
            w = min(tile.shape[1], plane.shape[1] - x0)
            plane[y0:y0 + h, x0:x0 + w] = tile[:h, :w]

    if true_size:
        scale = 2 ** level
        w = max(1, -(-int(true_size[0]) // scale))   # ceil
        h = max(1, -(-int(true_size[1]) // scale))
        plane = plane[:h, :w]

    return plane, level


def _signal(plane: np.ndarray, mask: np.ndarray, percentile: float) -> tuple[float, float]:
    """Return ``(signal_above_background, background)`` for one channel."""
    finite = plane[np.isfinite(plane)]
    if finite.size == 0:
        return float("nan"), float("nan")
    background = float(np.percentile(finite, BACKGROUND_PERCENTILE))
    inside = plane[mask & np.isfinite(plane)]
    if inside.size == 0:
        return float("nan"), background
    return float(np.percentile(inside, percentile)) - background, background


def measure_section(
    ets_path: str,
    exposures: Mapping[str, float],
    *,
    level: int | None = None,
    pixel_size_um: float | None = None,
    true_size: tuple[int, int] | None = None,
) -> SectionStain:
    """Measure the Aß and GFAP indices for one section.

    Args:
        ets_path: the tissue ``.ets`` container for this section.
        exposures: ``{"DAPI": ms, "FITC": ms, "Cy3": ms}``. Required -- exposure
            varies up to 12.6x here, so an unnormalised index is not comparable
            between animals.
        true_size: ``(width_px, height_px)`` at base level. Pass it -- see
            :func:`read_pyramid_plane` for why measuring on the padded grid
            inflates every signal.
    """
    out = SectionStain(ets_path=ets_path)
    if true_size is None:
        out.problems.append(
            "true_size not given - measuring on the padded tile grid would inflate "
            "the signal; refusing rather than returning a plausible wrong number"
        )
        return out
    try:
        dapi, lv = read_pyramid_plane(
            ets_path, CH_DAPI, level=level, pixel_size_um=pixel_size_um, true_size=true_size
        )
        cy3, _ = read_pyramid_plane(ets_path, CH_CY3, level=lv, true_size=true_size)
        fitc, _ = read_pyramid_plane(ets_path, CH_FITC, level=lv, true_size=true_size)
    except (ValueError, OSError) as exc:
        out.problems.append(str(exc))
        return out

    out.level = lv
    out.plane_width_px = dapi.shape[1]

    # Tissue from the counterstain, never from a marker channel: thresholding on Cy3
    # would define tissue as "where the thing we are measuring is", which guarantees
    # the answer.
    finite_dapi = dapi[np.isfinite(dapi)]
    if finite_dapi.size == 0:
        out.problems.append("DAPI plane is entirely missing tiles")
        return out
    # Otsu, but on the histogram clipped at its 99.5th percentile. Unclipped it is
    # not robust here: a handful of ultra-bright debris pixels (one section peaks at
    # 37,612 counts against a normal maximum near 670) pull the threshold up to
    # 8,889, so the "tissue" it finds is the dust speck -- 0.03% of the section --
    # and every statistic taken inside that mask describes an artefact. Clipping
    # costs nothing, since the threshold belongs far below the tissue peak anyway.
    ceiling = float(np.percentile(finite_dapi, 99.5))
    try:
        from skimage.filters import threshold_otsu
        thresh = float(threshold_otsu(np.minimum(finite_dapi, ceiling)))
    except Exception:  # pragma: no cover - skimage is a hard dependency
        thresh = float(np.percentile(finite_dapi, 75))
    mask = np.isfinite(dapi) & (dapi > thresh)
    out.tissue_fraction = float(mask.sum() / max(np.isfinite(dapi).sum(), 1))

    if out.tissue_fraction < 0.005:
        out.problems.append(f"almost no tissue found ({out.tissue_fraction:.2%})")
        return out

    dapi_sig, _ = _signal(dapi, mask, 50.0)
    cy3_sig, _ = _signal(cy3, mask, CY3_PERCENTILE)
    fitc_sig, _ = _signal(fitc, mask, CY3_PERCENTILE)

    e_dapi = float(exposures.get("DAPI") or float("nan"))
    e_cy3 = float(exposures.get("Cy3") or float("nan"))
    e_fitc = float(exposures.get("FITC") or float("nan"))

    denom = dapi_sig / e_dapi if dapi_sig and e_dapi else float("nan")
    if not denom or not np.isfinite(denom):
        out.problems.append("DAPI signal is zero or unmeasurable - cannot normalise")
        return out

    out.abeta_index = (cy3_sig / e_cy3) / denom
    out.gfap_index = (fitc_sig / e_fitc) / denom
    return out


def audit_manifest(df: Any, *, level: int | None = None, progress: Any = None) -> Any:
    """Measure every section in the manifest that has pixels, and compare.

    Returns a DataFrame with one row per measured section, carrying the recorded
    ``condition`` alongside the measured ``abeta_index``. Whether they agree is
    the caller's to judge and a human's to settle -- see the module docstring.
    """
    import pandas as pd

    rows = []
    for _, r in df.iterrows():
        path = r.get("ets_path")
        if not isinstance(path, str) or not path or not os.path.exists(path):
            continue
        w, h = r.get("width_px"), r.get("height_px")
        m = measure_section(
            path,
            {
                "DAPI": r.get("exposure_DAPI_ms"),
                "FITC": r.get("exposure_FITC_ms"),
                "Cy3": r.get("exposure_Cy3_ms"),
            },
            level=level,
            pixel_size_um=r.get("pixel_size_um"),
            true_size=(int(w), int(h)) if pd.notna(w) and pd.notna(h) else None,
        )
        rows.append({
            "tube_id": r.get("tube_id"),
            "scan": r.get("scan"),
            "section_label": r.get("section_label"),
            "condition": r.get("condition"),
            "use_for_measurement": r.get("use_for_measurement"),
            "abeta_index": m.abeta_index,
            "gfap_index": m.gfap_index,
            "tissue_fraction": m.tissue_fraction,
            "level": m.level,
            "plane_width_px": m.plane_width_px,
            "problems": "; ".join(m.problems),
        })
        if progress is not None:
            progress(rows[-1])
    return pd.DataFrame(rows)
