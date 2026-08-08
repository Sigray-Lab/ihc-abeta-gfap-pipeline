#!/usr/bin/env python3
"""The section manifest: one row per imaged brain section, records joined to metadata.

Why this module exists
----------------------
Everything downstream -- blinding, region delineation, classifier application,
aggregation -- reads this table and nothing else. It is the single place where the
*wet-lab record* (``config/slides.csv``) meets the *scanner metadata* (the ``.vsi``
index), and the single place where a section acquires the one property that decides
whether its pixels are data or a negative control.

The join in one paragraph
-------------------------
A slide holds three or four sections from one animal, arranged in two PAP-pen boxes.
The two sections inside a box share a staining condition: one box received primary
antibody, the other received DAPI + secondary only. **Box membership** is geometry and
comes from the stage coordinates (:func:`ihc.ingest.vsi_meta.assign_boxes`); the section
number is acquisition order and says nothing about slide position -- tube 49's positive
sections are ``01`` and ``04``. **Which box is positive** is a bench fact and comes only
from ``config/slides.csv``. Neither half can be recovered from the other, and neither
may be inferred from the images.

The three states of ``condition``
--------------------------------
``positive``
    The section sits in the box named by ``slides.csv:positive_box``, or that column
    says ``both`` (the slide was double-stained to salvage material and therefore has
    **no negative control**).
``negative``
    The section sits in the other box: DAPI + secondary only.
``unresolved``
    The ``slides.csv`` row carries a non-empty ``needs_confirmation``, or box
    membership could not be derived, or the slide has no ``slides.csv`` row at all.
    Every section of such a slide is excluded from the analysis manifest
    (``analysis_include == False``) and reported loudly. It is never guessed.
    Tube 37 is the live case: its ``positive_box`` column and its bench annotation
    name opposite boxes.

What a row is
-------------
Normally one imaged tissue series: ``row_kind == "section"``. Two bookkeeping kinds
exist so that a table read back from CSV can still be validated on its own:

``row_kind == "slides_csv_only"``
    A tube in ``slides.csv`` with no ``.vsi`` index file. One row, no section.
``in_slides_csv == False`` on a section row
    A ``.vsi`` with no ``slides.csv`` row -- condition is ``unresolved``.

Rescans
-------
``RawData/Rescan/`` holds re-acquisitions of tubes 51 and 60 (ADR-0006), two tissue
series each, **the positive box only**. The rescan renumbers its sections ``01``/``02``
regardless of which sections they were on the original slide, so joining a rescan's
``section_label`` to the original's would bind the wrong physical section and, on a
``far_label`` slide, the wrong condition. Rescan sections are therefore matched back to
the original scan **by stage X**, nearest-neighbour, with an explicit tolerance
(:data:`RESCAN_MATCH_TOLERANCE_FRACTION` of the original slide's smallest inter-section
spacing). Observed offsets in this cohort are 82-349 um against a smallest spacing of
6.7 mm, so the match is unambiguous by more than an order of magnitude -- but it is
checked rather than assumed, and a failed match yields ``condition == "unresolved"``.

Two consequences that anything consuming this table must handle:

*``(tube_id, section_label)`` is not a key.* Tubes 51 and 60 each carry two scans, and
the rescan's ``01``/``02`` collide with the original's ``01``/``02`` while naming
different physical sections with, on these two ``far_label`` slides, the **opposite**
condition. The key is ``(tube_id, scan, section_label)``. ``physical_section_label``
gives the section its identity on the original slide, so ``(tube_id,
physical_section_label)`` pairs the two scans of one section -- which is what the
imaging-repeatability check in ADR-0012 needs.

*Pick one scan before aggregating.* ``scan_is_preferred`` marks it, following
``config.yaml:intensity.exposure_correction.prefer_rescan``. Anything that collapses
the manifest to one row per section -- a blinded manifest, a QuPath project, an
animal-level mean -- must filter on it, or tube 51 contributes its low-exposure
original and its rescan as if they were two independent sections.

Missing pixels are normal
-------------------------
Payload folders have been transferred for 8 of 31 animals. "Metadata known, pixels not
yet on disk" is an ordinary state, not an error: ``payload_present`` is False, the row
is still built, and ``analysis_include`` is unaffected. Only pixel work filters on it.

Public API
----------
:func:`build_manifest`, :func:`write_manifest`, :func:`validate_manifest`,
:func:`crosscheck_condition_against_pixels`, :func:`format_slide_summary`.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import Any, Sequence

import numpy as np
import pandas as pd

from ihc.ingest.verify import (
    PIXEL_TYPES,
    find_companion,
    read_ets_index,
)
from ihc.ingest.vsi_meta import (
    MARGINAL_GAP_RATIO,
    BoxAssignmentError,
    VsiParseError,
    assign_boxes,
    read_vsi_meta,
)
from ihc.util.config import load_config, load_paths

__all__ = [
    "SCHEMA",
    "TOOL_VERSION",
    "MANIFEST_COLUMNS",
    "ManifestError",
    "build_manifest",
    "read_manifest",
    "write_manifest",
    "validate_manifest",
    "crosscheck_condition_against_pixels",
    "format_slide_summary",
    "format_crosscheck",
]

SCHEMA = "ihc.ingest.manifest/1"
TOOL_VERSION = "1.0.0"

# --------------------------------------------------------------------------
# Cohort expectations. Everything tunable is here so a reader can see in one
# place what "normal" is assumed to be, and so nothing is buried in a branch.
# --------------------------------------------------------------------------

#: Exposure the cohort was acquired at, milliseconds per channel. 29 of 31 slides
#: match; tubes 51 and 60 deviate and were rescanned (ADR-0006). Compared with a
#: tolerance because the scanner writes microseconds and we divide by 1000.
STANDARD_EXPOSURE_MS: dict[str, float] = {
    "DAPI": 128.547,
    "FITC": 397.927,
    "Cy3": 1839.999,
}
EXPOSURE_TOL_MS = 0.5

#: Channels as acquired, in order. The ETS channel axis is indexed in this order.
CHANNELS: tuple[str, ...] = ("DAPI", "FITC", "Cy3")

#: Treatment group (verbatim from ``slides.csv``) to delivery arm. The two arms are
#: analysed as two pre-specified comparisons, never pooled (spec section 10), so the
#: manifest has to carry the distinction explicitly rather than leaving it to a
#: substring test in every downstream script.
GROUP_TO_ARM: dict[str, str] = {
    "Rapamycin Diet": "diet",
    "Extra Control Diet": "diet",
    "Control IP (vehicle)": "ip",
    "Rapamycin IP": "ip",
}

#: A rescan section is matched to an original section by stage X. The match must be
#: closer than this fraction of the original slide's smallest inter-section spacing.
#: Observed worst case in this cohort is 349 um against 6722 um, i.e. 0.052.
RESCAN_MATCH_TOLERANCE_FRACTION = 0.25

#: Pixel size expected from the 20x objective, and the tolerance for a warning.
#: Never compare for bit-equality: the scanner writes a different value per series.
EXPECTED_PIXEL_SIZE_UM = 0.325
PIXEL_SIZE_TOL_UM = 0.001

#: Stack IDs that are never tissue.
NON_TISSUE_STACK_IDS = frozenset({1, 10000})

#: Column order of the manifest. Fixed, because a stable column order is what makes a
#: CSV diff readable in review.
MANIFEST_COLUMNS: tuple[str, ...] = (
    # identity
    "row_kind",
    "tube_id",
    "group",
    "arm",
    "scan",
    "section_label",
    "physical_section_label",
    "series_name",
    # the load-bearing join
    "box",
    "condition",
    "positive_box",
    "has_negative_control",
    "n_positive_sections_on_slide",
    "analysis_include",
    # slide-level record
    "slide_number_depth_index",
    "n_sections_on_slide",
    "n_sections_recorded",
    "needs_confirmation",
    "annotation",
    # geometry
    "stage_x_um",
    "stage_y_um",
    "box_gap_ratio",
    "box_split_gap_mm",
    "pixel_size_um",
    "width_px",
    "height_px",
    "area_px",
    "area_um2",
    # acquisition
    "exposure_DAPI_ms",
    "exposure_FITC_ms",
    "exposure_Cy3_ms",
    "exposure_is_standard",
    "acquisition_time",
    # provenance / storage
    "in_slides_csv",
    "in_vsi_index",
    "payload_present",
    "scan_is_preferred",
    "stack_id",
    "ets_path",
    "ets_sha256",
    "vsi_path",
    "vsi_sha256",
    "tube_id_in_file",
    "rescan_matched_section_label",
    "rescan_match_offset_um",
    # free text
    "notes",
    "warnings",
)

#: Nullable dtypes, applied after construction so that a missing integer stays missing
#: instead of turning the whole column into float (a section label silently becoming
#: 1.0 is exactly the class of bug the spec asks the schema to prevent).
_DTYPES: dict[str, str] = {
    "tube_id": "Int64",
    "slide_number_depth_index": "Int64",
    "n_sections_on_slide": "Int64",
    "n_sections_recorded": "Int64",
    "n_positive_sections_on_slide": "Int64",
    "width_px": "Int64",
    "height_px": "Int64",
    "area_px": "Int64",
    "stack_id": "Int64",
    "exposure_is_standard": "boolean",
    "has_negative_control": "boolean",
    "analysis_include": "boolean",
    "in_slides_csv": "boolean",
    "in_vsi_index": "boolean",
    "payload_present": "boolean",
    "scan_is_preferred": "boolean",
}

_STRING_COLUMNS = tuple(
    c for c in MANIFEST_COLUMNS
    if c not in _DTYPES
    and c not in {
        "stage_x_um", "stage_y_um", "box_gap_ratio", "box_split_gap_mm",
        "pixel_size_um", "area_um2", "exposure_DAPI_ms", "exposure_FITC_ms",
        "exposure_Cy3_ms", "rescan_match_offset_um",
    }
)

_FILENAME_TUBE_RE = re.compile(r"(\d+)")

#: A section label is ``"01"``..``"04"`` -- a two-character *string*, never a number.
#: ``pandas.read_csv`` will happily turn it into the integer 1, which then fails to
#: join against anything and cannot be told apart from a real value. This is the exact
#: failure the spec asks the schema to prevent, so it is checked rather than assumed.
_SECTION_LABEL_RE = re.compile(r"^0[1-9]$")


class ManifestError(RuntimeError):
    """The manifest could not be built at all (missing roster, unreadable config)."""


# ==========================================================================
# slides.csv -- the wet-lab record
# ==========================================================================


def _read_slides_csv(path: str | os.PathLike) -> dict[int, dict[str, Any]]:
    """Read ``config/slides.csv`` into ``{tube_id: row}``.

    Only structural problems raise here (a missing file, a duplicate tube). Semantic
    problems -- a ``positive_box`` that contradicts the geometry, a section count that
    disagrees with the index -- are the business of :func:`validate_manifest`, because
    the manifest must still be *built* so that the problem can be shown in context.
    """
    path = str(path)
    if not os.path.isfile(path):
        raise ManifestError(
            f"slides.csv not found at {path}. This file is the sole authority for which "
            f"PAP-pen box received primary antibody; without it no section has a "
            f"condition and no manifest can be built."
        )
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    rows: dict[int, dict[str, Any]] = {}
    for record in frame.to_dict("records"):
        raw = str(record.get("tube_id", "")).strip()
        if not raw:
            continue
        try:
            tube = int(raw)
        except ValueError as exc:
            raise ManifestError(f"slides.csv: tube_id {raw!r} is not an integer") from exc
        if tube in rows:
            raise ManifestError(f"slides.csv: duplicate tube_id {tube}")
        rows[tube] = {k: (v.strip() if isinstance(v, str) else v) for k, v in record.items()}
    if not rows:
        raise ManifestError(f"slides.csv at {path} has no usable rows")
    return rows


def _arm_for_group(group: str) -> str:
    """Map a treatment group to its delivery arm, or ``""`` if unrecognised."""
    return GROUP_TO_ARM.get(group, "")


def _int_or_none(value: Any) -> int | None:
    try:
        text = str(value).strip()
        return int(text) if text else None
    except (TypeError, ValueError):
        return None


# ==========================================================================
# payload discovery
# ==========================================================================


def _tissue_stacks(vsi_path: str) -> list[tuple[int, str]]:
    """Return ``[(stack_id, ets_path), ...]`` for the tissue stacks of one slide.

    Ascending stack ID. Empty when the payload folder is absent, which is the normal
    state for 23 of the 31 animals.
    """
    companion = find_companion(vsi_path)
    if companion is None:
        return []
    out: list[tuple[int, str]] = []
    try:
        entries = sorted(os.listdir(companion))
    except OSError:
        return []
    for entry in entries:
        if not (entry.startswith("stack") and entry[5:].isdigit()):
            continue
        stack_id = int(entry[5:])
        if stack_id in NON_TISSUE_STACK_IDS:
            continue
        ets = os.path.join(companion, entry, "frame_t_0.ets")
        if os.path.isfile(ets):
            out.append((stack_id, ets))
    return sorted(out)


def _bind_sections_to_stacks(
    section_labels: Sequence[str],
    stacks: Sequence[tuple[int, str]],
    true_sizes: dict[str, tuple[int | None, int | None]],
) -> tuple[dict[str, tuple[int, str]], list[str]]:
    """Bind each section label to its stack folder, and check the binding.

    Olympus writes tissue stacks in acquisition order (``stack10002``, ``10005``,
    ``10008``, ``10011``) and the section number *is* acquisition order, so ascending
    stack ID pairs with ascending section label. That is a convention, and binding the
    wrong pixels to a section is the worst failure this pipeline can have -- the run
    completes and every number lands in a plausible range -- so the pairing is
    confirmed against the true image rectangle from VSI tag 2053: the tile grid must be
    able to hold the section and must not exceed it by more than one tile in each
    direction.

    Returns:
        ``({section_label: (stack_id, ets_path)}, problems)``. On any inconsistency the
        binding is returned **empty** and the problem is reported: no pixels are better
        than the wrong pixels.
    """
    problems: list[str] = []
    labels = sorted(section_labels)
    if not stacks:
        return {}, problems
    if len(stacks) != len(labels):
        problems.append(
            f"the index names {len(labels)} tissue series but the payload holds "
            f"{len(stacks)} tissue stack(s) {[s for s, _ in stacks]}; not binding "
            f"pixels to sections for this slide"
        )
        return {}, problems

    binding: dict[str, tuple[int, str]] = {}
    for label, (stack_id, ets) in zip(labels, stacks):
        binding[label] = (stack_id, ets)
        width, height = true_sizes.get(label, (None, None))
        if width is None or height is None:
            continue
        try:
            index = read_ets_index(ets)
        except Exception as exc:  # noqa: BLE001 - reported, never fatal
            problems.append(f"stack{stack_id}: could not read the tile index ({exc})")
            continue
        base_level = min(c[0][-1] for c in index["chunks"]) if index["chunks"] else None
        if base_level is None:
            problems.append(f"stack{stack_id}: tile index is empty")
            continue
        base = [c for c in index["chunks"] if c[0][-1] == base_level]
        cols = max(c[0][0] for c in base) + 1
        rows = max(c[0][1] for c in base) + 1
        tile_w, tile_h, _ = index["tile"]
        # The grid may fall short of the true image (whole edge tile columns are absent
        # from the scanner sample mask -- tube 30 does this) but it may never need more
        # than one extra tile beyond it.
        if cols * tile_w > width + tile_w or rows * tile_h > height + tile_h:
            problems.append(
                f"stack{stack_id} bound to section _{label}: its {cols}x{rows} grid of "
                f"{tile_w}x{tile_h} tiles cannot belong to a {width}x{height} px section. "
                f"The payload folder may belong to a different slide. Refusing to bind "
                f"pixels for this slide."
            )
            return {}, problems
    return binding, problems


# ==========================================================================
# verification hashes
# ==========================================================================


def _verified_hashes(paths: dict[str, Any] | None = None) -> dict[str, dict[str, Any]]:
    """Read SHA-256 digests recorded by ``./ihc verify``, keyed by absolute ``.vsi`` path.

    ``vsi_sha256`` and ``ets_sha256`` are null in a freshly built manifest and are
    filled in from the verification reports when they exist, so the manifest can state
    *whether* the bytes it describes have been checked rather than implying they have.
    Absent or unreadable reports are not an error.
    """
    if paths is None:
        try:
            paths = load_paths()
        except Exception:  # noqa: BLE001 - config problems are reported elsewhere
            return {}
    out: dict[str, dict[str, Any]] = {}
    verify_dir = paths.get("work.verify_dir")
    if verify_dir is None:
        return out
    for name in ("verify_raw.json", "verify_rescan.json"):
        report = os.path.join(str(verify_dir), name)
        if not os.path.isfile(report):
            continue
        try:
            with open(report, encoding="utf-8") as handle:
                summary = json.load(handle)
        except (OSError, ValueError):
            continue
        for dataset in summary.get("datasets", []):
            vsi = dataset.get("vsi_path")
            if not vsi:
                continue
            per_section = {
                str(stack.get("section")): stack.get("sha256")
                for stack in dataset.get("stacks", [])
                if stack.get("role") == "tissue" and stack.get("section")
            }
            out[os.path.abspath(vsi)] = {
                "vsi_sha256": dataset.get("vsi_sha256"),
                "ets_sha256": per_section,
                "verified_ok": bool(dataset.get("ok")),
            }
    return out


# ==========================================================================
# rescan matching
# ==========================================================================


def _side_of_box_split(x_um, original_series):
    """Which PAP-pen box does stage position *x_um* fall in?

    Only meaningful because the geometry is so lopsided: the gap BETWEEN the two
    boxes is 11-23 mm in this cohort while sections within a box sit 6-8 mm apart.
    A position is assigned only if it sits clearly on one side of that gap and not
    inside it; otherwise None, and the caller refuses rather than guessing.
    """
    xs = sorted(s.stage_x_um for s in original_series)
    if len(xs) < 2:
        return None
    gaps = [(b - a, a, b) for a, b in zip(xs, xs[1:])]
    widest, lo, hi = max(gaps)
    others = [g for g, _, _ in gaps if (g, lo, hi) != (widest, lo, hi)]
    # The split must actually be a split: clearly wider than the within-box spacing.
    if others and widest < 1.5 * max(others):
        return None
    if lo < x_um < hi:
        return None                       # inside the gap itself -- no safe call
    return "near_label" if x_um <= lo else "far_label"


def _match_rescan_to_original(
    rescan_series: Sequence[Any], original_series: Sequence[Any]
) -> tuple[dict[str, tuple[str, float]], list[str]]:
    """Match rescan sections to original sections by stage X, nearest neighbour.

    The rescan renumbers its sections ``01``/``02`` whatever they were on the original
    slide, so the label cannot be joined directly. Stage X can: the slide goes back on
    the same stage in the same orientation, and the sections are millimetres apart.

    Every match must be closer than :data:`RESCAN_MATCH_TOLERANCE_FRACTION` of the
    original slide's smallest inter-section spacing, and no two rescan sections may
    claim the same original. A failure returns an empty mapping rather than a guess.

    Returns:
        ``({rescan_label: (original_label, offset_um)}, problems)``
    """
    problems: list[str] = []
    if not rescan_series or len(original_series) < 2:
        return {}, ["cannot match a rescan without at least two original sections"]

    xs = sorted(s.stage_x_um for s in original_series)
    spacing = min(b - a for a, b in zip(xs, xs[1:]))
    tolerance = RESCAN_MATCH_TOLERANCE_FRACTION * spacing

    mapping: dict[str, tuple[str, float]] = {}
    claimed: dict[str, str] = {}
    newly_scanned: dict[str, str] = {}
    for section in rescan_series:
        best = min(original_series, key=lambda o: abs(o.stage_x_um - section.stage_x_um))
        offset = abs(best.stage_x_um - section.stage_x_um)
        if offset > tolerance:
            # No original at this position. Two readings, and they are distinguishable.
            #
            # Either the slide moved between sessions (in which case EVERY section is
            # offset and nothing matches), or the second session imaged a section the
            # first one skipped. The wet lab confirmed the latter for tubes 33, 42 and
            # 54: "because I had some extra time on the booking slot I put it also the
            # ones I had one positive". Those slides always had four sections; the
            # first scan only captured three.
            #
            # A genuinely new section still has an unambiguous box, because the gap
            # BETWEEN boxes is 11-23 mm while sections within a box sit 6-8 mm apart.
            # So place it by which side of the box split it falls on, and mark it
            # newly_scanned so nothing downstream mistakes it for a re-acquisition.
            box_side = _side_of_box_split(section.stage_x_um, original_series)
            if box_side is not None:
                newly_scanned[section.section_label] = box_side
                problems.append(
                    f"rescan section _{section.section_label} is {offset / 1000:.2f} mm "
                    f"from any original section, which is beyond the "
                    f"{tolerance / 1000:.2f} mm matching tolerance. Read as a section "
                    f"the first scan did not capture, not a re-acquisition; placed in "
                    f"the {box_side} box by stage position. (Confirmed at the bench for "
                    f"tubes 33, 42 and 54.)"
                )
                continue
            problems.append(
                f"rescan section _{section.section_label} is {offset / 1000:.2f} mm from "
                f"the nearest original section (_{best.section_label}); the tolerance is "
                f"{tolerance / 1000:.2f} mm (a quarter of the {spacing / 1000:.2f} mm "
                f"smallest inter-section spacing), and it does not fall cleanly on "
                f"either side of the box split. Condition cannot be carried over."
            )
            return {}, problems
        if best.section_label in claimed:
            problems.append(
                f"rescan sections _{claimed[best.section_label]} and "
                f"_{section.section_label} both match original section "
                f"_{best.section_label}; the match is ambiguous"
            )
            return {}, problems
        claimed[best.section_label] = section.section_label
        mapping[section.section_label] = (best.section_label, offset)
    return mapping, problems, newly_scanned


# ==========================================================================
# build
# ==========================================================================


def _blank_row() -> dict[str, Any]:
    return {column: pd.NA for column in MANIFEST_COLUMNS}


def _exposure_is_standard(exposure_ms: dict[str, float]) -> bool:
    """True if all three channels match the cohort's standard exposure."""
    for channel, expected in STANDARD_EXPOSURE_MS.items():
        value = exposure_ms.get(channel)
        if value is None or abs(value - expected) > EXPOSURE_TOL_MS:
            return False
    return True


def _slide_rows(
    vsi_path: str,
    scan: str,
    slides: dict[int, dict[str, Any]],
    hashes: dict[str, dict[str, Any]],
    original_meta: dict[int, Any],
    prefer_rescan: bool,
    has_rescan: set[int],
) -> list[dict[str, Any]]:
    """Build every manifest row for one ``.vsi`` index file."""
    rows: list[dict[str, Any]] = []
    notes: list[str] = []

    try:
        meta = read_vsi_meta(vsi_path)
    except VsiParseError as exc:
        row = _blank_row()
        tube_match = _FILENAME_TUBE_RE.search(os.path.basename(vsi_path))
        row.update(
            row_kind="section",
            tube_id=int(tube_match.group(1)) if tube_match else pd.NA,
            scan=scan,
            condition="unresolved",
            analysis_include=False,
            in_slides_csv=False,
            in_vsi_index=False,
            payload_present=False,
            vsi_path=os.path.abspath(vsi_path),
            notes="the .vsi index could not be read",
            warnings=str(exc),
        )
        return [row]

    tube = meta.tube_id
    record = slides.get(tube, {}) if tube is not None else {}
    in_slides_csv = bool(record)
    if not in_slides_csv:
        notes.append(
            "no config/slides.csv row for this tube, so the positive box is unknown"
        )

    positive_box = record.get("positive_box", "")
    needs_confirmation = record.get("needs_confirmation", "")
    # Treatment group comes from the CUSTODIAN tree, never from the committed
    # slides.csv -- see load_group_allocation(). Anyone without the custodian file
    # gets an empty group and a blinded manifest, which is the correct default for
    # everyone who is not the custodian.
    group = _GROUPS.get(tube, record.get("group", "")) if tube is not None else ""
    arm = _arm_for_group(group)
    if group and not arm:
        notes.append(f"group {group!r} is not one of {sorted(GROUP_TO_ARM)}")

    # --- box membership, from stage X and nothing else --------------------
    boxes: dict[str, Any] | None = None
    box_of: dict[str, str] = {}
    rescan_match: dict[str, tuple[str, float]] = {}
    if scan == "rescan":
        original = original_meta.get(tube)
        if original is None:
            notes.append(
                "no original scan for this rescan, so box membership cannot be "
                "carried over"
            )
        else:
            rescan_match, match_problems, newly_scanned = _match_rescan_to_original(
                meta.series, original.series
            )
            notes.extend(match_problems)
            if rescan_match or newly_scanned:
                try:
                    boxes = assign_boxes(original.series)
                except BoxAssignmentError as exc:
                    notes.append(f"box assignment failed on the original scan: {exc}")
                else:
                    for label, (original_label, _offset) in rescan_match.items():
                        box_of[label] = (
                            "near_label"
                            if original_label in boxes["near_label"]
                            else "far_label"
                        )
                    # Sections the first scan never captured: box comes from stage
                    # position, not from a matched original. Confirmed at the bench
                    # for tubes 33, 42 and 54.
                    for label, side in newly_scanned.items():
                        box_of[label] = side
                        # It has no counterpart in the first scan, so it needs its own
                        # physical identity. Reusing the rescan's own label would
                        # collide with the original section that happens to share it
                        # (tube 33: rescan _01 is a new section, but original _01 is a
                        # different piece of tissue), and the two would be merged.
                        rescan_match[label] = (f"new{label}", 0.0)
    else:
        try:
            boxes = assign_boxes(meta.series)
        except BoxAssignmentError as exc:
            notes.append(f"box assignment failed: {exc}")
        else:
            for label in boxes["near_label"]:
                box_of[label] = "near_label"
            for label in boxes["far_label"]:
                box_of[label] = "far_label"
            if boxes["gap_ratio"] < MARGINAL_GAP_RATIO:
                notes.append(
                    f"between-box stage gap is only {boxes['gap_ratio']:.2f}x the "
                    f"largest within-box gap; the split is real but tight"
                )

    # --- condition ---------------------------------------------------------
    # THE LOAD-BEARING DERIVATION. Records first, geometry second, pixels never.
    unresolved_reason = ""
    if not in_slides_csv:
        unresolved_reason = "no slides.csv row"
    elif needs_confirmation:
        unresolved_reason = f"slides.csv needs_confirmation: {needs_confirmation}"
    elif positive_box not in ("near_label", "far_label", "both"):
        unresolved_reason = f"slides.csv positive_box is {positive_box!r}"
    elif positive_box != "both" and not box_of:
        unresolved_reason = "box membership could not be derived from stage X"

    def condition_for(label: str) -> str:
        if unresolved_reason:
            return "unresolved"
        if positive_box == "both":
            return "positive"
        box = box_of.get(label)
        if box is None:
            return "unresolved"
        return "positive" if box == positive_box else "negative"

    conditions = {s.section_label: condition_for(s.section_label) for s in meta.series}
    n_positive = sum(1 for v in conditions.values() if v == "positive")

    # --- pixels on disk ----------------------------------------------------
    stacks = _tissue_stacks(vsi_path)
    binding, bind_problems = _bind_sections_to_stacks(
        [s.section_label for s in meta.series],
        stacks,
        {s.section_label: (s.width_px, s.height_px) for s in meta.series},
    )
    notes.extend(bind_problems)

    digest = hashes.get(os.path.abspath(vsi_path), {})
    ets_digests = digest.get("ets_sha256", {})

    slide_warnings = list(meta.warnings)
    if tube is not None and tube in has_rescan:
        preferred = (scan == "rescan") if prefer_rescan else (scan == "original")
    else:
        preferred = True

    for series in meta.series:
        label = series.section_label
        row = _blank_row()
        stack_id, ets = binding.get(label, (None, None))
        matched_label, matched_offset = rescan_match.get(label, ("", None))
        width, height = series.width_px, series.height_px
        area_px = width * height if (width and height) else None
        pixel_size = series.pixel_size_um
        row.update(
            row_kind="section",
            tube_id=tube if tube is not None else pd.NA,
            group=group,
            arm=arm,
            scan=scan,
            section_label=label,
            physical_section_label=matched_label or label,
            series_name=series.name,
            box=box_of.get(label, ""),
            condition=conditions[label],
            positive_box=positive_box,
            # filled in animal-wide by build_manifest once every scan has been read
            has_negative_control=pd.NA,
            n_positive_sections_on_slide=n_positive if not unresolved_reason else pd.NA,
            analysis_include=conditions[label] != "unresolved",
            slide_number_depth_index=_int_or_none(record.get("slide_number_depth_index")),
            n_sections_on_slide=meta.n_tissue_series,
            n_sections_recorded=_int_or_none(record.get("n_sections")),
            needs_confirmation=needs_confirmation,
            annotation=record.get("annotation", ""),
            stage_x_um=round(series.stage_x_um, 3),
            stage_y_um=round(series.stage_y_um, 3),
            box_gap_ratio=round(boxes["gap_ratio"], 4) if boxes else None,
            box_split_gap_mm=round(boxes["split_gap_mm"], 4) if boxes else None,
            pixel_size_um=pixel_size,
            width_px=width,
            height_px=height,
            area_px=area_px,
            area_um2=(
                round(area_px * pixel_size * pixel_size, 1)
                if area_px and pixel_size
                else None
            ),
            exposure_DAPI_ms=series.exposure_ms.get("DAPI"),
            exposure_FITC_ms=series.exposure_ms.get("FITC"),
            exposure_Cy3_ms=series.exposure_ms.get("Cy3"),
            exposure_is_standard=_exposure_is_standard(series.exposure_ms),
            acquisition_time=series.acquisition_time,
            in_slides_csv=in_slides_csv,
            in_vsi_index=True,
            payload_present=ets is not None,
            scan_is_preferred=preferred,
            stack_id=stack_id if stack_id is not None else pd.NA,
            ets_path=ets or "",
            ets_sha256=ets_digests.get(label) or "",
            vsi_path=os.path.abspath(vsi_path),
            vsi_sha256=digest.get("vsi_sha256") or "",
            tube_id_in_file=meta.tube_id_in_file or "",
            rescan_matched_section_label=matched_label,
            rescan_match_offset_um=(
                round(matched_offset, 1) if matched_offset is not None else None
            ),
            notes="; ".join(
                n for n in ([unresolved_reason] if unresolved_reason else []) + notes
            ),
            warnings="; ".join(slide_warnings),
        )
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# Per-section notes from the bench
# --------------------------------------------------------------------------- #

def load_section_notes(path=None):
    """Read ``config/section_notes.csv``: per-section facts only the bench knows.

    slides.csv answers one question per *slide* — which PAP-pen box got antibody.
    This file carries the finer-grained things that turn up when someone actually
    looks at the images: a hole in one section, dirt on another, which of two scans
    is the better one. They matter because they change what gets measured, and they
    are not recoverable from the pixels alone (a hole and a fold look similar; dirt
    and a plaque both look bright).

    Columns: tube_id, scan, section_label, action, note, source, date.
    ``section_label`` may be blank when the action applies to a whole scan.
    Recognised actions:
      exclude_section  -- do not measure this section at all
      flag_artefact    -- measurable, but something must be excluded by the artefact mask
      prefer_section   -- the better member of a positive pair
      prefer_scan      -- use this scan of this tube, overriding the cohort default

    Returns a list of dicts. A missing file is normal and returns [].
    """
    import csv as _csv
    if path is None:
        path = _repo_root() / "config" / "section_notes.csv"
    path = Path(path)
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return [
            {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            for row in _csv.DictReader(fh)
            if row.get("tube_id", "").strip()
        ]


def _repo_root():
    return Path(__file__).resolve().parents[3]


def apply_section_notes(df, notes=None):
    """Fold ``config/section_notes.csv`` into a built manifest.

    Adds three columns and may flip one:
      bench_note          -- the free text, so it travels with the row
      bench_action        -- exclude_section / flag_artefact / prefer_section / ""
      has_flagged_artefact-- True where the bench says something must be masked out
      scan_is_preferred   -- OVERRIDDEN by a prefer_scan note

    The prefer_scan override matters. The cohort default prefers a rescan, which is
    right for tubes 51 and 60 — they were re-acquired specifically to fix exposure.
    It is wrong for tube 49, where the rescan's section 01 was out of focus and the
    original is the better image. A global rule cannot know that; the bench can.
    """
    pd_ = _pandas() if "_pandas" in globals() else __import__("pandas")
    if notes is None:
        notes = load_section_notes()
    df = df.copy()
    for col, default in (("bench_note", ""), ("bench_action", "")):
        if col not in df.columns:
            df[col] = default
    df["has_flagged_artefact"] = False
    df["flagged_lower_quality"] = False

    for note in notes:
        tube = int(note["tube_id"])
        scan = (note.get("scan") or "").strip()
        label = (note.get("section_label") or "").strip()
        action = (note.get("action") or "").strip()

        if action == "prefer_scan":
            # whole-tube override: this scan preferred, every other scan of it not
            same_tube = df["tube_id"] == tube
            if not same_tube.any():
                continue
            df.loc[same_tube, "scan_is_preferred"] = (df.loc[same_tube, "scan"] == scan)
            df.loc[same_tube, "bench_note"] = note.get("note", "")
            # Mark the action too, so finalise_scan_preference can honour the override.
            # Without this the note was invisible to the later recomputation and tube 49
            # silently reverted to the rescan the bench had rejected as out of focus.
            df.loc[same_tube & (df["scan"] == scan), "bench_action"] = "prefer_scan"
            continue

        sel = (df["tube_id"] == tube)
        if scan:
            sel &= (df["scan"] == scan)
        if label:
            sel &= (df["section_label"].astype(str) == label)
        if not sel.any():
            continue
        df.loc[sel, "bench_note"] = note.get("note", "")
        df.loc[sel, "bench_action"] = action
        if action == "flag_artefact":
            df.loc[sel, "has_flagged_artefact"] = True
        elif action == "flag_quality":
            # Usable but the bench rates it lower. Kept in the analysis; carried as a
            # flag so a sensitivity analysis without these sections is one filter away.
            if "flagged_lower_quality" not in df.columns:
                df["flagged_lower_quality"] = False
            df.loc[sel, "flagged_lower_quality"] = True
        elif action == "exclude_section":
            if "analysis_include" in df.columns:
                df.loc[sel, "analysis_include"] = False
    return df


def finalise_scan_preference(df):
    """Recompute ``scan_is_preferred`` now that conditions are known, and derive the
    ONE column later stages should filter on.

    This exists because of a real defect. ``scan_is_preferred`` was computed from the
    mere *existence* of a file in ``Rescan/``, before anything knew whether that rescan
    had resolved. For tubes 33, 42 and 54 the "rescan" turned out to contain a section
    the original does not have, so the pair failed to match, the rescan rows came out
    unresolved -- and the originals had already been demoted. The animal then had no
    usable row at all, and no validator noticed.

    Two filters existed and both were wrong: ``analysis_include`` alone double-counted
    tubes 49, 51 and 60 (importing the very exposures the rescans exist to replace),
    while ``scan_is_preferred`` alone deleted tubes 33, 42 and 54 outright -- three
    animals from three different groups. Both wrong answers looked plausible.

    Rule: a rescan is preferred only if it actually yielded usable rows for that tube.
    Otherwise the original stands. ``use_for_measurement`` is then the single column
    every later stage filters on.
    """
    df = df.copy()
    usable = df["condition"].isin(["positive", "negative"]) & df["analysis_include"]

    for tube, grp in df.groupby("tube_id"):
        # A bench note can pin the whole tube to one scan (tube 49: the rescan was out
        # of focus, so the original stands even though a rescan exists).
        pinned = None
        if "bench_action" in grp:
            ov = grp.loc[grp["bench_action"] == "prefer_scan", "scan"]
            if len(ov):
                pinned = ov.iloc[0]

        for _sec, sgrp in grp.groupby("physical_section_label", dropna=False):
            good = sgrp[usable[sgrp.index]]
            if pinned is not None and (sgrp["scan"] == pinned).any():
                winner_idx = sgrp.index[sgrp["scan"] == pinned]
            elif (good["scan"] == "rescan").any():
                winner_idx = good.index[good["scan"] == "rescan"]
            elif len(good):
                winner_idx = good.index[good["scan"] == sorted(set(good["scan"]))[0]]
            else:
                winner_idx = sgrp.index[:1]
            df.loc[sgrp.index, "scan_is_preferred"] = False
            df.loc[winner_idx[:1], "scan_is_preferred"] = True

    df["use_for_measurement"] = df["analysis_include"] & df["scan_is_preferred"] & usable
    return df


def load_group_allocation(path=None):
    """Read tube -> treatment group from the CUSTODIAN tree, not from the repo.

    The allocation table used to live in ``config/slides.csv``, which is committed.
    That put tube-to-group in plain text in every clone, fork and reflog of the
    repository -- an unblinding shortcut that needed no key at all, and one that
    ``.gitignore`` was carefully protecting the key from while publishing the same
    information as CSV. It now lives beside the blinding key.

    Returns {} if the file is absent, which is the correct state for anyone who is
    not the custodian: the manifest then carries an empty ``group`` and the pipeline
    runs blinded, as it should.
    """
    import csv as _csv
    if path is None:
        try:
            from ihc.util.config import load_paths
            path = load_paths()["custodian_root"] / "group_allocation.csv"
        except Exception:
            return {}
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {int(r["tube_id"]): r["group"] for r in _csv.DictReader(fh)
                if r.get("tube_id", "").strip()}


_GROUPS: dict = {}


def inputs_fingerprint(raw_root=None, slides_csv=None) -> dict:
    """A cheap fingerprint of everything the manifest is built from.

    Exists because derived artefacts went stale twice without anyone noticing. Payload
    folders arrive in batches over weeks, and a QuPath project built last Tuesday looks
    exactly like one built this morning -- it just quietly contains fewer animals. The
    second time, the delineation project was short by seven animals and nothing said so.

    Deliberately coarse: the set of index files, which of them have pixels, and a digest
    of the two records that decide condition. That is enough to detect "the inputs moved"
    without hashing 43 GB.
    """
    import hashlib
    from pathlib import Path as _P
    if raw_root is None:
        from ihc.util.config import load_paths
        paths = load_paths()
        raw_root = paths["raw_root"]
        slides_csv = slides_csv or paths["config.slides_csv"]
    raw_root = _P(raw_root)

    vsi = sorted(p.name for p in raw_root.glob("*.vsi"))
    with_pixels = sorted(
        p.name for p in raw_root.glob("_Image_*_") if any(p.rglob("*.ets")))
    rescan_dir = raw_root / "Rescan"
    rescans = sorted(p.name for p in rescan_dir.glob("*.vsi")) if rescan_dir.is_dir() else []

    h = hashlib.sha256()
    for item in (vsi, with_pixels, rescans):
        h.update("\n".join(item).encode())
    for record in (slides_csv, _P(slides_csv).parent / "section_notes.csv" if slides_csv else None):
        if record and _P(record).exists():
            h.update(_P(record).read_bytes())

    return {
        "n_vsi": len(vsi),
        "n_payloads_with_pixels": len(with_pixels),
        "n_rescans": len(rescans),
        "digest": h.hexdigest()[:16],
    }


def build_manifest(
    raw_root: str | os.PathLike | None = None,
    slides_csv: str | os.PathLike | None = None,
    *,
    include_rescans: bool = True,
) -> pd.DataFrame:
    """Build the section manifest: one row per imaged tissue section.

    Joins the scanner metadata in every ``.vsi`` index under *raw_root* to the wet-lab
    staining record in *slides_csv*, deriving PAP-pen box membership from the stage
    coordinates and ``condition`` from the record. No pixels are read.

    Args:
        raw_root: Folder holding the ``.vsi`` index files. Defaults to
            ``paths.yaml:roots.raw_root``. Subfolders are not searched except for the
            rescan folder, which is handled explicitly.
        slides_csv: The wet-lab record. Defaults to ``paths.yaml:config.slides_csv``.
        include_rescans: Also read ``RawData/Rescan/``. The rescans are re-acquisitions
            of tubes 51 and 60 at standard exposure (ADR-0006) and carry only the
            positive box, so their sections are matched back to the original scan by
            stage X before a condition is assigned.

    Returns:
        A :class:`pandas.DataFrame` with the columns of :data:`MANIFEST_COLUMNS`,
        sorted by ``(tube_id, scan, section_label)``. Rows whose ``condition`` is
        ``unresolved`` are present but carry ``analysis_include == False``.

    Raises:
        ManifestError: *slides_csv* is missing or structurally unusable, or *raw_root*
            does not exist. An individual unreadable ``.vsi`` never raises: it becomes
            a row with ``in_vsi_index == False`` and a note.

    Example:
        >>> df = build_manifest()                                    # doctest: +SKIP
        >>> df.groupby("condition").size()                           # doctest: +SKIP
        condition
        negative       50
        positive       68
        unresolved      4
    """
    global _GROUPS
    _GROUPS = load_group_allocation()
    paths = load_paths()
    root = str(raw_root) if raw_root is not None else str(paths["raw_root"])
    csv_path = str(slides_csv) if slides_csv is not None else str(paths["config.slides_csv"])
    if not os.path.isdir(root):
        raise ManifestError(f"raw_root does not exist: {root}")

    slides = _read_slides_csv(csv_path)
    hashes = _verified_hashes(paths)
    try:
        prefer_rescan = bool(
            load_config()["intensity"]["exposure_correction"]["prefer_rescan"]
        )
    except (KeyError, TypeError):
        prefer_rescan = True

    vsi_glob = str(paths["raw"]["vsi_glob"])
    originals = sorted(
        os.path.join(root, name)
        for name in os.listdir(root)
        if name.lower().endswith(vsi_glob.lstrip("*"))
    )

    rescan_dir = os.path.join(root, str(paths["raw"]["rescan_subdir"]))
    rescans: list[str] = []
    if include_rescans and os.path.isdir(rescan_dir):
        rescans = sorted(
            os.path.join(rescan_dir, name)
            for name in os.listdir(rescan_dir)
            if name.lower().endswith(vsi_glob.lstrip("*"))
        )

    # The rescans need their original scan's geometry to recover box membership, so
    # the originals are parsed first and kept.
    original_meta: dict[int, Any] = {}
    for path in originals:
        try:
            meta = read_vsi_meta(path)
        except VsiParseError:
            continue
        if meta.tube_id is not None:
            original_meta[meta.tube_id] = meta

    has_rescan: set[int] = set()
    for path in rescans:
        match = _FILENAME_TUBE_RE.search(os.path.basename(path))
        if match:
            has_rescan.add(int(match.group(1)))

    rows: list[dict[str, Any]] = []
    for path in originals:
        rows.extend(
            _slide_rows(path, "original", slides, hashes, original_meta,
                        prefer_rescan, has_rescan)
        )
    for path in rescans:
        rows.extend(
            _slide_rows(path, "rescan", slides, hashes, original_meta,
                        prefer_rescan, has_rescan)
        )

    # A tube in slides.csv with no .vsi at all gets one bookkeeping row, so that the
    # roster survives a CSV round-trip and validate_manifest can still see it.
    seen = {row["tube_id"] for row in rows if row["tube_id"] is not pd.NA}
    for tube in sorted(set(slides) - seen):
        record = slides[tube]
        row = _blank_row()
        row.update(
            row_kind="slides_csv_only",
            tube_id=tube,
            group=_GROUPS.get(tube, record.get("group", "")),
            arm=_arm_for_group(_GROUPS.get(tube, record.get("group", ""))),
            scan="original",
            section_label="",
            condition="unresolved",
            positive_box=record.get("positive_box", ""),
            analysis_include=False,
            slide_number_depth_index=_int_or_none(record.get("slide_number_depth_index")),
            n_sections_recorded=_int_or_none(record.get("n_sections")),
            needs_confirmation=record.get("needs_confirmation", ""),
            annotation=record.get("annotation", ""),
            in_slides_csv=True,
            in_vsi_index=False,
            payload_present=False,
            scan_is_preferred=False,
            notes="slides.csv names this tube but there is no .vsi index file for it",
        )
        rows.append(row)

    frame = pd.DataFrame(rows, columns=list(MANIFEST_COLUMNS))

    # `has_negative_control` is a property of the ANIMAL, not of one scan. The rescans
    # carry only the positive box, so a per-scan flag would say "no negative control"
    # for tube 51 while two negative sections of tube 51 sit on the original slide.
    # The negative-control gate is defined per animal-where-negatives-exist (D-10), so
    # animal-level is the semantics every consumer of this column wants.
    for tube, group in frame.groupby("tube_id", dropna=True):
        sections = group[group["row_kind"] == "section"]
        if sections.empty or (sections["condition"] == "unresolved").all():
            continue
        frame.loc[group.index, "has_negative_control"] = bool(
            (sections["condition"] == "negative").any()
        )

    for column, dtype in _DTYPES.items():
        frame[column] = frame[column].astype(dtype)
    for column in _STRING_COLUMNS:
        frame[column] = frame[column].fillna("").astype(str)
    frame = frame.sort_values(
        ["tube_id", "scan", "section_label"], kind="stable"
    ).reset_index(drop=True)
    frame.attrs["schema"] = SCHEMA
    frame.attrs["tool_version"] = TOOL_VERSION
    frame.attrs["built_utc"] = _dt.datetime.now(_dt.timezone.utc).isoformat(
        timespec="seconds"
    )
    frame.attrs["raw_root"] = os.path.abspath(root)
    frame.attrs["slides_csv"] = os.path.abspath(csv_path)
    return frame


# ==========================================================================
# reading back
# ==========================================================================


def read_manifest(path: str | os.PathLike) -> pd.DataFrame:
    """Read ``manifest.csv`` back with the dtypes it was written with.

    Use this rather than a bare :func:`pandas.read_csv`. Left to itself pandas infers
    ``section_label`` ``"01"`` as the integer ``1`` and ``stack_id`` ``10002`` as the
    float ``10002.0``. The section label then joins against nothing, and a label that
    has become a number is indistinguishable from a label that was always wrong --
    which is the specific silent failure the spec asks the manifest schema to prevent.
    """
    dtypes = {column: "string" for column in _STRING_COLUMNS}
    dtypes.update(
        {column: dtype for column, dtype in _DTYPES.items() if dtype != "boolean"}
    )
    frame = pd.read_csv(path, dtype=dtypes, keep_default_na=True)
    for column, dtype in _DTYPES.items():
        if dtype == "boolean" and column in frame:
            frame[column] = frame[column].astype("boolean")
    for column in _STRING_COLUMNS:
        if column in frame:
            frame[column] = frame[column].fillna("").astype(str)
    return frame


# ==========================================================================
# validation
# ==========================================================================


def validate_manifest(df: pd.DataFrame) -> list[str]:
    """Return every problem found in *df*, as sentences a human can act on.

    An empty list means the manifest is internally consistent **and** consistent with
    the wet-lab record. It does not mean the science is right: that is what the
    negative-control gate and the pixel auditor are for.

    Checked here, in the order the list is returned:

    1. a section whose ``condition`` is ``unresolved`` (tube 37 today);
    2. a slide whose section count disagrees between ``slides.csv`` and the VSI index;
    3. a slide with no positive section;
    4. a tube in ``slides.csv`` with no ``.vsi``, or a ``.vsi`` with no ``slides.csv`` row;
    5. a tube whose filename and in-file ID (VSI tag 2061) disagree;
    6. a duplicate ``(tube_id, scan, section_label)``;
    7. supporting checks: box assignment failures, ``has_negative_control`` disagreeing
       with the section counts, non-standard exposure, pixel size out of tolerance,
       tube 59 present, an unrecognised treatment group, and a slide whose pixels could
       not be bound to its sections.

    Args:
        df: A manifest from :func:`build_manifest`, or one read back from
            ``manifest.csv``. Both work: every check reads columns only.

    Returns:
        Problems, most structural first. Ordinary states are **not** problems:
        a missing payload folder, a three-section slide, a slide with no negative
        control because both boxes were stained.
    """
    problems: list[str] = []
    if df.empty:
        return ["the manifest is empty"]

    missing = [column for column in MANIFEST_COLUMNS if column not in df.columns]
    if missing:
        return [f"the manifest is missing required column(s): {missing}"]

    sections = df[df["row_kind"] == "section"].copy()

    # 0. dtype integrity. A section label is a two-character string; if it has arrived
    # as an integer the file was read with inferred dtypes and every join downstream
    # will fail quietly.
    bad_labels = [
        str(v) for v in sections["section_label"].dropna().unique()
        if str(v) and not _SECTION_LABEL_RE.match(str(v))
    ]
    if bad_labels:
        problems.append(
            f"section_label holds values that are not zero-padded two-character strings: "
            f"{sorted(bad_labels)[:8]}. This normally means the CSV was read with "
            f"inferred dtypes ('01' becomes the integer 1). Read it with "
            f"ihc.ingest.manifest.read_manifest() instead."
        )
        # Carry on so the remaining checks still run and the report is complete.
        sections["section_label"] = sections["section_label"].astype(str)

    # 1. unresolved --------------------------------------------------------
    # Sections from an unreadable index are unresolved too, but they get their own
    # message below; listing them twice buries the ones a human can actually act on.
    unresolved = sections[
        (sections["condition"] == "unresolved") & sections["in_vsi_index"].fillna(False)
    ]
    for tube, group in unresolved.groupby("tube_id", dropna=False):
        reason = next((n for n in group["notes"] if n), "no reason recorded")
        problems.append(
            f"tube {tube}: condition is UNRESOLVED for {len(group)} section(s) "
            f"({', '.join(sorted(group['section_label']))}) -- {reason}. These sections "
            f"are excluded from the analysis manifest and must not be guessed; the "
            f"answer comes from the bench, not from the images."
        )

    # 2. section count disagreement ---------------------------------------
    originals = sections[sections["scan"] == "original"]
    for tube, group in originals.groupby("tube_id", dropna=False):
        recorded = group["n_sections_recorded"].dropna().unique()
        if len(recorded) == 0:
            continue
        recorded_n = int(recorded[0])
        indexed_n = int(group["n_sections_on_slide"].dropna().iloc[0])
        if recorded_n != indexed_n:
            problems.append(
                f"tube {tube}: slides.csv records {recorded_n} section(s) but the VSI "
                f"index holds {indexed_n}. One of the two is wrong -- a truncated index "
                f"loses whole sections silently, and a wrong count in the record means "
                f"the box split is being read against the wrong layout."
            )

    # 3. no positive section ----------------------------------------------
    for (tube, scan), group in sections.groupby(["tube_id", "scan"], dropna=False):
        if (group["condition"] == "unresolved").all():
            continue  # already reported above
        if not (group["condition"] == "positive").any():
            problems.append(
                f"tube {tube} ({scan}): no positive section. Every slide must contribute "
                f"at least one antibody-stained section; a slide of nothing but negative "
                f"controls contributes no data and is almost certainly a wrong "
                f"positive_box in slides.csv."
            )

    # 4. roster mismatches -------------------------------------------------
    for tube in sorted(df.loc[df["row_kind"] == "slides_csv_only", "tube_id"].dropna()):
        problems.append(
            f"tube {tube}: named in slides.csv but there is no .vsi index file for it. "
            f"Either the file was never transferred or the record names an animal that "
            f"was not imaged."
        )
    unreadable = sections[~sections["in_vsi_index"].fillna(False)]
    for _, row in unreadable.iterrows():
        problems.append(
            f"{os.path.basename(str(row['vsi_path'])) or 'a .vsi'}: the index could not "
            f"be read ({row['warnings'] or 'no detail'}). A truncated index loses whole "
            f"sections silently and yields a confidently wrong box assignment; "
            f"re-transfer the file."
        )
    orphan = sections[
        ~sections["in_slides_csv"].fillna(False) & sections["in_vsi_index"].fillna(False)
    ]
    for tube in sorted(orphan["tube_id"].dropna().unique()):
        problems.append(
            f"tube {tube}: a .vsi index exists but slides.csv has no row for it, so the "
            f"positive box is unknown and every section is unresolved."
        )

    # 5. filename vs in-file tube ID --------------------------------------
    # One message per FILE, not per section: the mismatch is a property of the file and
    # four identical lines make the list harder to read, not more convincing.
    for path, group in sections.groupby("vsi_path", dropna=False):
        row = group.iloc[0]
        in_file = str(row.get("tube_id_in_file") or "").strip()
        if not in_file or pd.isna(row["tube_id"]):
            continue
        try:
            matches = int(in_file) == int(row["tube_id"])
        except ValueError:
            matches = False
        if not matches:
            problems.append(
                f"{os.path.basename(str(path))}: the filename says tube "
                f"{row['tube_id']} but VSI tag 2061 says {in_file!r}. The file may have "
                f"been renamed, or the wrong slide was scanned under this name. Which of "
                f"the two is right decides whose brain these pixels are."
            )

    # 6. duplicates --------------------------------------------------------
    key = ["tube_id", "scan", "section_label"]
    duplicated = sections[sections.duplicated(key, keep=False)]
    for values, group in duplicated.groupby(key, dropna=False):
        problems.append(
            f"duplicate rows for tube {values[0]} {values[1]} section _{values[2]} "
            f"({len(group)} rows). Every section must appear exactly once."
        )

    # 7. supporting checks -------------------------------------------------
    for (tube, scan), group in sections.groupby(["tube_id", "scan"], dropna=False):
        if (group["box"] == "").all() and (group["condition"] != "unresolved").any():
            problems.append(
                f"tube {tube} ({scan}): no PAP-pen box could be assigned to any section, "
                f"yet a condition was still derived. That should be impossible."
            )

    # has_negative_control is an ANIMAL-level flag: the rescans carry only the positive
    # box, so it is checked across every scan of a tube rather than within one.
    for tube, group in sections.groupby("tube_id", dropna=False):
        flags = group["has_negative_control"].dropna().unique()
        if len(flags) > 1:
            problems.append(
                f"tube {tube}: has_negative_control takes more than one value "
                f"({sorted(bool(f) for f in flags)}) across its rows; it is a property "
                f"of the animal and must be constant."
            )
            continue
        if len(flags) != 1:
            continue
        expected = bool((group["condition"] == "negative").any())
        if bool(flags[0]) != expected:
            problems.append(
                f"tube {tube}: has_negative_control={bool(flags[0])} but the manifest "
                f"holds {int((group['condition'] == 'negative').sum())} negative "
                f"section(s) for this animal."
            )

    # Non-standard exposure is a problem only when it is *unmitigated*. Tubes 51 and 60
    # were acquired off-standard and rescanned (ADR-0006), so their originals are a
    # pre-specified sensitivity analysis rather than an open defect; listing a handled
    # fact as a problem every run is how a problem list stops being read.
    non_standard = sections[sections["exposure_is_standard"] == False]  # noqa: E712
    rescanned = set(sections.loc[sections["scan"] == "rescan", "tube_id"].dropna())
    for tube in sorted(non_standard["tube_id"].dropna().unique()):
        rows_for_tube = non_standard[non_standard["tube_id"] == tube]
        scans = sorted(set(rows_for_tube["scan"]))
        if tube in rescanned and scans == ["original"]:
            continue
        problems.append(
            f"tube {tube} ({', '.join(scans)}): acquired at NON-STANDARD exposure and "
            f"NOT rescanned. Dividing by exposure rescales but does not restore photons "
            f"-- on tube 51 the tissue sits ~7 grey levels above background where on "
            f"tube 29 it sits ~208 -- so a numerical correction is a sensitivity "
            f"analysis, not a fix (ADR-0006). Rescan this slide."
        )

    off_pixel = sections[
        sections["pixel_size_um"].notna()
        & ((sections["pixel_size_um"] - EXPECTED_PIXEL_SIZE_UM).abs() > PIXEL_SIZE_TOL_UM)
    ]
    for _, row in off_pixel.iterrows():
        problems.append(
            f"tube {row['tube_id']} section _{row['section_label']}: pixel size "
            f"{row['pixel_size_um']:.6f} um is more than {PIXEL_SIZE_TOL_UM} um from the "
            f"expected {EXPECTED_PIXEL_SIZE_UM} um -- wrong objective, or the wrong slide."
        )

    if 59 in set(df["tube_id"].dropna().astype(int)):
        problems.append(
            "tube 59 is present. It was excluded before imaging (mounting fault, PAP-pen "
            "and antibody leak) and must not appear anywhere in the manifest."
        )

    unknown_arm = sections[(sections["group"] != "") & (sections["arm"] == "")]
    for tube in sorted(unknown_arm["tube_id"].dropna().unique()):
        group_name = unknown_arm.loc[unknown_arm["tube_id"] == tube, "group"].iloc[0]
        problems.append(
            f"tube {tube}: treatment group {group_name!r} does not map to an arm; the "
            f"two arms are analysed as two pre-specified comparisons and cannot be pooled."
        )

    for _, row in sections.iterrows():
        note = str(row.get("notes") or "")
        if "Refusing to bind pixels" in note or "not binding" in note:
            problems.append(
                f"tube {row['tube_id']} ({row['scan']}): {note}"
            )
            break

    tube_col = "tube_id"
    # --- R1 invariant: every tube must survive the filter later stages use ----------
    # Both previous filters silently lost or duplicated animals and nothing complained.
    if "use_for_measurement" in df.columns:
        for tube, grp in df.groupby(tube_col):
            keep = grp[grp["use_for_measurement"]]
            if keep.empty:
                problems.append(
                    f"tube {tube}: NO row survives use_for_measurement -- this animal "
                    f"would silently vanish from the analysis. Scans present: "
                    f"{sorted(set(grp['scan']))}.")
                continue
            if not (keep["condition"] == "positive").any():
                problems.append(
                    f"tube {tube}: rows survive the filter but none is a positive "
                    f"section, so the animal contributes no measurable data.")
            dupes = keep.duplicated(subset=[tube_col, "physical_section_label"]).sum()
            if dupes:
                problems.append(
                    f"tube {tube}: {dupes} physical section(s) appear more than once "
                    f"after filtering -- that animal would be double-counted.")

    return problems


# ==========================================================================
# writing
# ==========================================================================


def write_manifest(df: pd.DataFrame, out_dir: str | os.PathLike) -> dict[str, str]:
    """Write the manifest and its analysis subset to *out_dir*.

    Four files, all regenerable:

    ``manifest.csv`` / ``manifest.json``
        Every row, including the unresolved ones. This is the provenance record.
    ``manifest_analysis.csv``
        Only ``analysis_include == True`` -- the table later stages read. Slides whose
        condition is unresolved are absent from it by construction, so a stage cannot
        accidentally quantify a section nobody has resolved.
    ``manifest_problems.txt``
        The output of :func:`validate_manifest`, so the problems travel with the table.

    Args:
        df: A manifest from :func:`build_manifest`.
        out_dir: Destination, created if absent. Normally ``work_root/manifest``.

    Returns:
        ``{"manifest_csv": ..., "manifest_json": ..., "analysis_csv": ...,
        "problems_txt": ...}`` as absolute paths.
    """
    out = os.path.abspath(str(out_dir))
    os.makedirs(out, exist_ok=True)

    manifest_csv = os.path.join(out, "manifest.csv")
    manifest_json = os.path.join(out, "manifest.json")
    analysis_csv = os.path.join(out, "manifest_analysis.csv")
    problems_txt = os.path.join(out, "manifest_problems.txt")

    df.to_csv(manifest_csv, index=False)

    payload = {
        "schema": df.attrs.get("schema", SCHEMA),
        "tool_version": df.attrs.get("tool_version", TOOL_VERSION),
        "built_utc": df.attrs.get(
            "built_utc",
            _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        ),
        "raw_root": df.attrs.get("raw_root", ""),
        "slides_csv": df.attrs.get("slides_csv", ""),
        "n_rows": int(len(df)),
        "n_sections": int((df["row_kind"] == "section").sum()),
        "n_analysis_rows": int(df["analysis_include"].fillna(False).sum()),
        "counts_by_condition": {
            str(k): int(v)
            for k, v in df.loc[df["row_kind"] == "section", "condition"]
            .value_counts()
            .items()
        },
        "columns": list(MANIFEST_COLUMNS),
        "rows": json.loads(df.to_json(orient="records")),
    }
    with open(manifest_json, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    analysis = df[df["analysis_include"].fillna(False)]
    analysis.to_csv(analysis_csv, index=False)

    problems = validate_manifest(df)
    with open(problems_txt, "w", encoding="utf-8") as handle:
        if problems:
            handle.write(f"{len(problems)} problem(s)\n\n")
            for problem in problems:
                handle.write(f"- {problem}\n")
        else:
            handle.write("no problems\n")

    return {
        "manifest_csv": manifest_csv,
        "manifest_json": manifest_json,
        "analysis_csv": analysis_csv,
        "problems_txt": problems_txt,
    }


# ==========================================================================
# human-readable summary
# ==========================================================================


def format_slide_summary(df: pd.DataFrame) -> str:
    """Render one line per slide: layout, condition split, exposure, pixels."""
    sections = df[df["row_kind"] == "section"]
    header = (
        f"{'tube':>4} {'scan':<8} {'group':<20} {'arm':<4} {'n':>1} "
        f"{'pos box':<10} {'positive':<12} {'negative':<12} {'negctl':<6} "
        f"{'exp':<4} {'pixels':<7} {'condition':<10}"
    )
    lines = [header, "-" * len(header)]
    for (tube, scan), group in sections.groupby(["tube_id", "scan"], dropna=False):
        tube = "?" if pd.isna(tube) else int(tube)
        group = group.sort_values("section_label")
        positive = ",".join(group.loc[group["condition"] == "positive", "section_label"])
        negative = ",".join(group.loc[group["condition"] == "negative", "section_label"])
        unresolved = (group["condition"] == "unresolved").any()
        flag = group["has_negative_control"].dropna()
        neg_ctl = "-" if flag.empty else ("yes" if bool(flag.iloc[0]) else "NO")
        exposure = "std" if bool(group["exposure_is_standard"].all()) else "NON"
        pixels = "yes" if bool(group["payload_present"].all()) else (
            "part" if bool(group["payload_present"].any()) else "index"
        )
        lines.append(
            f"{tube:>4} {scan:<8} {str(group['group'].iloc[0])[:20]:<20} "
            f"{str(group['arm'].iloc[0]):<4} {len(group):>1} "
            f"{str(group['positive_box'].iloc[0]):<10} {positive or '-':<12} "
            f"{negative or '-':<12} {neg_ctl:<6} {exposure:<4} {pixels:<7} "
            f"{'UNRESOLVED' if unresolved else 'ok':<10}"
        )
    for _, row in df[df["row_kind"] == "slides_csv_only"].iterrows():
        lines.append(
            f"{int(row['tube_id']):>4} {'-':<8} {str(row['group'])[:20]:<20} "
            f"{str(row['arm']):<4} 0 {str(row['positive_box']):<10} "
            f"{'-':<12} {'-':<12} {'-':<6} {'-':<4} {'none':<7} {'NO .vsi':<10}"
        )
    return "\n".join(lines)


# ==========================================================================
# ==========================================================================
#   AUDITOR -- image-based cross-check of the record. NEVER a source of truth.
# ==========================================================================
# ==========================================================================

#: Below this the two conditions are not separable and the section is reported
#: ``inconclusive`` rather than being forced into a call. The measured cohort values
#: are 0.57-1.41 for negative sections and 4.86-20.07 for positive ones across the 34
#: sections with pixels, so the band between 2.0 and 3.0 is empty by a wide margin.
GFAP_INDEX_NEGATIVE_MAX = 2.0
GFAP_INDEX_POSITIVE_MIN = 3.0

#: Fewer tissue pixels than this at the chosen pyramid level and the percentiles are
#: not trustworthy; the section is reported ``inconclusive``.
MIN_TISSUE_PX = 20_000

#: Pyramid level used by default. Level 3 is downsample 8 (2.6 um/px), which reads a
#: whole section in well under a second and leaves the 99.9th percentile intact.
DEFAULT_CROSSCHECK_LEVEL = 3

_CROSSCHECK_COLUMNS: tuple[str, ...] = (
    "tube_id",
    "scan",
    "section_label",
    "condition_from_record",
    "pixel_condition",
    "agreement",
    "gfap_index",
    "abeta_index",
    "within_slide_separation",
    "fitc_p999_tissue",
    "cy3_p999_tissue",
    "dapi_p99_support",
    "dapi_tissue_threshold",
    "tissue_fraction_of_support",
    "n_tissue_px",
    "n_support_px",
    "level",
    "ets_path",
    "note",
)


def _decode_tile(buffer: bytes, compression: int, dtype: str) -> "np.ndarray | None":
    """Decode one ETS tile. Returns ``None`` if the codec is not available."""
    import imagecodecs

    if compression == 3:
        return imagecodecs.jpeg2k_decode(buffer)
    if compression == 2:
        return imagecodecs.jpeg8_decode(buffer)
    if compression in (0, 5):
        return np.frombuffer(buffer, dtype=np.dtype(dtype))
    return None


def _read_level(ets_path: str, level: int) -> tuple["np.ndarray | None", "np.ndarray | None", str]:
    """Assemble one pyramid level of an ``.ets`` into ``(channels, support, note)``.

    The tile grid is **sparse** -- 4-18 % of the positions inside the bounding box were
    never acquired, because the scanner used a sample mask. Absent positions are left
    zero in the mosaic and False in the returned ``support`` mask, so they can be
    excluded from every statistic rather than counted as dark tissue.

    Returns ``(None, None, note)`` when the level is absent or nothing decodes.
    """
    try:
        index = read_ets_index(ets_path)
    except Exception as exc:  # noqa: BLE001 - reported, never fatal
        return None, None, f"could not read the tile index: {exc}"

    chunks = [c for c in index["chunks"] if c[0][-1] == level]
    if not chunks:
        levels = sorted({c[0][-1] for c in index["chunks"]})
        return None, None, f"pyramid level {level} is absent (present: {levels})"

    tile_w, tile_h, _ = index["tile"]
    cols = max(c[0][0] for c in chunks) + 1
    rows = max(c[0][1] for c in chunks) + 1
    channel_axis = 2 if index["n_dim"] >= 4 else None
    channels = (
        sorted({c[0][channel_axis] for c in chunks}) if channel_axis is not None else [0]
    )
    if len(channels) != len(CHANNELS):
        return None, None, (
            f"{len(channels)} channel plane(s), expected {len(CHANNELS)} "
            f"({', '.join(CHANNELS)})"
        )

    dtype, _bytes = PIXEL_TYPES.get(index["pixel_type"], ("uint16", 2))
    mosaic = np.zeros((len(channels), rows * tile_h, cols * tile_w), dtype=np.uint16)
    support = np.zeros((rows * tile_h, cols * tile_w), dtype=bool)
    decoded = 0
    with open(ets_path, "rb") as handle:
        for coords, offset, nbytes in chunks:
            x, y = coords[0], coords[1]
            channel = coords[channel_axis] if channel_axis is not None else 0
            handle.seek(offset)
            buffer = handle.read(nbytes)
            try:
                tile = _decode_tile(buffer, index["compression"], dtype)
            except Exception:  # noqa: BLE001 - one bad tile must not stop the section
                continue
            if tile is None:
                return None, None, (
                    f"compression code {index['compression']} is not decodable here"
                )
            if tile.ndim == 1:
                if tile.size != tile_w * tile_h:
                    continue
                tile = tile.reshape(tile_h, tile_w)
            height, width = tile.shape[:2]
            plane = channels.index(channel)
            mosaic[plane, y * tile_h : y * tile_h + height, x * tile_w : x * tile_w + width] = tile
            support[y * tile_h : y * tile_h + height, x * tile_w : x * tile_w + width] = True
            decoded += 1

    if decoded == 0:
        return None, None, "no tile decoded at this level"
    return mosaic, support, ""


def _section_pixel_stats(ets_path: str, level: int) -> dict[str, Any]:
    """Compute the DAPI-normalised GFAP and Abeta indices for one section.

    Tissue is defined from DAPI alone -- the nuclear counterstain is present in every
    section regardless of staining condition, so using it cannot make a negative
    control look like tissue-free glass. The threshold is the geometric midpoint
    between the 25th and 99.9th percentile of DAPI over the acquired support, which is
    scale-free and therefore survives the 2.1x DAPI exposure difference on tube 51.

    Both indices are then divided by the 99th percentile of DAPI over the same support.
    That is the normalisation: it cancels the per-slide gain that would otherwise let
    an exposure difference masquerade as a staining difference. The 99th percentile is
    used rather than the 99.9th because DAPI carries occasional saturated debris that
    would swamp the extreme tail.
    """
    mosaic, support, note = _read_level(ets_path, level)
    if mosaic is None:
        return {"note": note}

    dapi = mosaic[0].astype(np.float64)
    fitc = mosaic[1].astype(np.float64)
    cy3 = mosaic[2].astype(np.float64)
    values = dapi[support]
    if values.size < MIN_TISSUE_PX:
        return {"note": f"only {values.size} acquired pixels at level {level}"}

    low = max(float(np.percentile(values, 25)), 1.0)
    high = max(float(np.percentile(values, 99.9)), 2.0)
    threshold = math.sqrt(low * high)
    tissue = support & (dapi > threshold)
    n_tissue = int(tissue.sum())

    result: dict[str, Any] = {
        "dapi_tissue_threshold": round(threshold, 2),
        "n_support_px": int(support.sum()),
        "n_tissue_px": n_tissue,
        "tissue_fraction_of_support": round(n_tissue / max(int(support.sum()), 1), 4),
        "note": "",
    }
    if n_tissue < MIN_TISSUE_PX:
        result["note"] = (
            f"only {n_tissue} tissue pixels at level {level}; too few for a percentile"
        )
        return result

    dapi_p99 = max(float(np.percentile(values, 99)), 1.0)
    fitc_p999 = float(np.percentile(fitc[tissue], 99.9))
    cy3_p999 = float(np.percentile(cy3[tissue], 99.9))
    result.update(
        dapi_p99_support=round(dapi_p99, 1),
        fitc_p999_tissue=round(fitc_p999, 1),
        cy3_p999_tissue=round(cy3_p999, 1),
        gfap_index=round(fitc_p999 / dapi_p99, 3),
        abeta_index=round(cy3_p999 / dapi_p99, 3),
    )
    return result


def crosscheck_condition_against_pixels(
    df: pd.DataFrame,
    *,
    level: int = DEFAULT_CROSSCHECK_LEVEL,
    max_slides: int | None = None,
) -> pd.DataFrame:
    """AUDITOR ONLY. Ask the pixels whether they agree with ``condition``, and report.

    **This function never writes to ``condition`` and nothing in this module gives it
    a way to.** It returns a separate report. The reason is not stylistic:

    *Deriving condition from pixels would silently convert a failed positive stain into
    a fake negative control.* A section that was stained with primary antibody but did
    not take would be relabelled ``negative``, disappear from the numerator, and then be
    used to certify that the assay has no background -- which is precisely the thing it
    is evidence against.

    *And it would make the negative-control QC gate circular.* That gate asks whether
    the frozen classifier calls near-zero positive area on no-primary tissue. If the
    no-primary sections had been chosen *because* they are dim, the gate would be
    testing the selection rule rather than the classifier, and would pass by
    construction however bad the classifier was.

    So: **records are truth; pixels audit.** A disagreement is a blocking flag for
    human review (ADR-0003, spec section 2) and is resolved at the bench, never by
    the code.

    How the call is made
    --------------------
    Per section, at a low pyramid level: tissue is thresholded on DAPI, and the 99.9th
    percentile of GFAP (FITC) within tissue is divided by the 99th percentile of DAPI
    over the acquired support. Sparse, never-acquired tile positions are excluded
    throughout -- they are missing support, not background (ADR-0010). The section is
    called ``positive`` above :data:`GFAP_INDEX_POSITIVE_MIN`, ``negative`` below
    :data:`GFAP_INDEX_NEGATIVE_MAX`, and ``inconclusive`` in between. Measured across
    the 34 sections with pixels: negatives 0.57-1.41, positives 4.86-20.07, and the
    smallest within-slide separation is 8.0x.

    Args:
        df: A manifest from :func:`build_manifest`. Only rows with
            ``payload_present == True`` are read; the rest are skipped silently,
            because 23 of 31 animals are index-only and that is a normal state.
        level: Pyramid level. 0 is full resolution; 3 (the default) is downsample 8.
        max_slides: Stop after this many slides. For a quick look; the full cohort of
            payloads takes about a minute.

    Returns:
        One row per audited section with ``condition_from_record``,
        ``pixel_condition``, ``agreement`` (``agree`` / ``CONTRADICTS`` /
        ``inconclusive`` / ``not_checked``), both indices, the raw components they were
        built from, and ``within_slide_separation`` -- the smallest positive index
        divided by the largest negative index on that slide, or NA where the slide
        carries only one condition.
    """
    sections = df[
        (df["row_kind"] == "section")
        & df["payload_present"].fillna(False)
        & (df["ets_path"].astype(str) != "")
    ].copy()
    if sections.empty:
        return pd.DataFrame(columns=list(_CROSSCHECK_COLUMNS))

    slides = list(dict.fromkeys(zip(sections["tube_id"], sections["scan"])))
    if max_slides is not None:
        slides = slides[:max_slides]
    wanted = set(slides)

    rows: list[dict[str, Any]] = []
    for _, section in sections.iterrows():
        if (section["tube_id"], section["scan"]) not in wanted:
            continue
        stats = _section_pixel_stats(str(section["ets_path"]), level)
        record_condition = str(section["condition"])
        gfap = stats.get("gfap_index")
        if gfap is None:
            pixel_condition = "inconclusive"
        elif gfap >= GFAP_INDEX_POSITIVE_MIN:
            pixel_condition = "positive"
        elif gfap <= GFAP_INDEX_NEGATIVE_MAX:
            pixel_condition = "negative"
        else:
            pixel_condition = "inconclusive"

        if record_condition == "unresolved":
            agreement = "not_checked"
        elif pixel_condition == "inconclusive":
            agreement = "inconclusive"
        elif pixel_condition == record_condition:
            agreement = "agree"
        else:
            agreement = "CONTRADICTS"

        row = {column: pd.NA for column in _CROSSCHECK_COLUMNS}
        row.update(
            tube_id=section["tube_id"],
            scan=section["scan"],
            section_label=section["section_label"],
            condition_from_record=record_condition,
            pixel_condition=pixel_condition,
            agreement=agreement,
            level=level,
            ets_path=section["ets_path"],
        )
        row.update({k: v for k, v in stats.items() if k in _CROSSCHECK_COLUMNS})
        rows.append(row)

    report = pd.DataFrame(rows, columns=list(_CROSSCHECK_COLUMNS))
    if report.empty:
        return report

    # Within-slide separation: the cleanest evidence available, because both conditions
    # on one slide share staining batch, exposure and section thickness.
    for (tube, scan), group in report.groupby(["tube_id", "scan"], dropna=False):
        positive = group.loc[group["condition_from_record"] == "positive", "gfap_index"].dropna()
        negative = group.loc[group["condition_from_record"] == "negative", "gfap_index"].dropna()
        if positive.empty or negative.empty:
            continue
        separation = float(positive.min()) / max(float(negative.max()), 1e-9)
        mask = (report["tube_id"] == tube) & (report["scan"] == scan)
        report.loc[mask, "within_slide_separation"] = round(separation, 2)

    report["note"] = report["note"].fillna("")
    return report.sort_values(["tube_id", "scan", "section_label"]).reset_index(drop=True)


def format_crosscheck(report: pd.DataFrame) -> str:
    """Render the auditor's report as a table, contradictions last and loud."""
    if report.empty:
        return "  no sections with pixels on disk; nothing to audit"
    header = (
        f"{'tube':>4} {'scan':<8} {'sec':<3} {'record':<10} {'pixels':<12} "
        f"{'gfap_idx':>9} {'abeta_idx':>10} {'sep':>7}  verdict"
    )
    lines = [header, "-" * len(header)]
    def number(value: Any, width: int, places: int = 2) -> str:
        return f"{'-':>{width}}" if pd.isna(value) else f"{value:{width}.{places}f}"

    for _, row in report.iterrows():
        lines.append(
            f"{row['tube_id']:>4} {row['scan']:<8} {row['section_label']:<3} "
            f"{row['condition_from_record']:<10} {row['pixel_condition']:<12} "
            f"{number(row['gfap_index'], 9)} {number(row['abeta_index'], 10)} "
            f"{number(row['within_slide_separation'], 7, 1)}  {row['agreement']}"
            + (f"   [{row['note']}]" if row["note"] else "")
        )
    contradictions = report[report["agreement"] == "CONTRADICTS"]
    lines.append("")
    if contradictions.empty:
        lines.append(
            f"  {int((report['agreement'] == 'agree').sum())} of {len(report)} audited "
            f"section(s) agree with slides.csv; 0 contradictions."
        )
    else:
        lines.append(
            f"  {len(contradictions)} CONTRADICTION(S) between the record and the pixels. "
            f"This BLOCKS. It is resolved at the bench, never by relabelling: "
            f"see ADR-0003."
        )
    return "\n".join(lines)


# ==========================================================================
# command line
# ==========================================================================


def _main(argv: Sequence[str]) -> int:  # pragma: no cover - thin wrapper
    """``python -m ihc.ingest.manifest [--crosscheck] [--level N]``."""
    import argparse

    parser = argparse.ArgumentParser(prog="python -m ihc.ingest.manifest")
    parser.add_argument("--out", default=None, help="output directory")
    parser.add_argument("--crosscheck", action="store_true", help="run the pixel auditor")
    parser.add_argument("--level", type=int, default=DEFAULT_CROSSCHECK_LEVEL)
    parser.add_argument("--max-slides", type=int, default=None)
    parser.add_argument("--no-rescans", action="store_true")
    args = parser.parse_args(list(argv))

    frame = build_manifest(include_rescans=not args.no_rescans)
    out = args.out or str(load_paths()["work.manifest_dir"])
    written = write_manifest(frame, out)
    print(format_slide_summary(frame))
    print()
    problems = validate_manifest(frame)
    for problem in problems:
        print(f"PROBLEM  {problem}")
    print(f"\n{len(frame)} row(s) -> {written['manifest_csv']}")
    if args.crosscheck:
        report = crosscheck_condition_against_pixels(
            frame, level=args.level, max_slides=args.max_slides
        )
        print()
        print(format_crosscheck(report))
    return 1 if problems else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main(sys.argv[1:]))
