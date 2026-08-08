#!/usr/bin/env python3
"""Validate config/slides.csv. Plain pandas and stdlib — no schema framework.

slides.csv is the authoritative record of which PAP-pen box carries primary antibody.
Getting it wrong means quantifying negative controls as data, and every downstream
number would look entirely plausible. So the checks are blunt and loud.

Run:  python3 config/validate_config.py [path/to/slides.csv]
Exit: 0 clean (warnings allowed), 1 on any error.
"""
import sys
from pathlib import Path

import pandas as pd

HEADER = ("tube_id", "group", "slide_number_depth_index", "positive_box", "n_sections",
          "no_negative_control", "annotation", "image_crosscheck", "needs_confirmation")
BOXES = {"near_label", "far_label", "both"}
THREE_SECTION_TUBES = {30, 33, 34, 42, 53, 54}   # verified from the payloads


def validate_slides(csv_path):
    """Return (errors, warnings) as two lists of strings."""
    errors, warnings = [], []
    df = pd.read_csv(csv_path, dtype={"tube_id": "Int64", "n_sections": "Int64"})

    if tuple(df.columns) != HEADER:
        errors.append(f"header mismatch\n  expected {HEADER}\n  found    {tuple(df.columns)}")
        return errors, warnings   # nothing below is meaningful with the wrong columns

    for row in df.itertuples(index=False):
        where = f"tube {row.tube_id}"
        if pd.isna(row.tube_id) or not (29 <= int(row.tube_id) <= 60):
            errors.append(f"{where}: tube_id outside 29..60")
        if row.positive_box not in BOXES:
            errors.append(f"{where}: positive_box {row.positive_box!r} not in {sorted(BOXES)}")
        if pd.isna(row.n_sections) or int(row.n_sections) not in (3, 4):
            errors.append(f"{where}: n_sections {row.n_sections} not in (3, 4)")
        elif int(row.tube_id) in THREE_SECTION_TUBES and int(row.n_sections) != 3:
            errors.append(f"{where}: known three-section animal recorded as {row.n_sections}")
        elif int(row.tube_id) not in THREE_SECTION_TUBES and int(row.n_sections) == 3:
            errors.append(f"{where}: recorded as 3 sections but is not a known three-section animal")
        if (row.positive_box == "both") != (str(row.no_negative_control).strip().lower() == "yes"):
            errors.append(f"{where}: positive_box={row.positive_box} contradicts "
                          f"no_negative_control={row.no_negative_control}")
        if not pd.isna(row.needs_confirmation) and str(row.needs_confirmation).strip():
            warnings.append(f"{where}: NEEDS CONFIRMATION -- {str(row.needs_confirmation).strip()}")

    dupes = df.loc[df.tube_id.duplicated(keep=False), "tube_id"].dropna().unique().tolist()
    if dupes:
        errors.append(f"duplicate tube_id: {dupes}")
    if 59 in set(df.tube_id.dropna().astype(int)):
        errors.append("tube 59 present: excluded before imaging (mounting fault), must not appear")
    if len(df) != 31:
        warnings.append(f"{len(df)} rows; the cohort is 31 slides (tubes 29-58 and 60)")
    return errors, warnings


def crosscheck_annotations_against_geometry(csv_path, raw_root):
    """Do the section numbers named in a row's free-text annotation actually fall
    inside the box that `positive_box` points at?

    This exists because `slides.csv` is the sole authority for which PAP-pen box got
    primary antibody, and a wrong `positive_box` means two negative-control sections
    get quantified as data while every number stays in a plausible range. The bench
    annotation is an independent statement of the same fact, so the two can be made
    to check each other -- and on the first run this caught tube 37, whose box column
    and annotation name opposite boxes.

    Costs ~0.15 s over all 31 index files. Returns a list of warning strings.
    """
    import csv as _csv
    import re as _re
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from ihc.ingest.vsi_meta import read_vsi_meta, assign_boxes

    out = []
    for row in _csv.DictReader(open(csv_path)):
        annotation = row.get("annotation") or ""
        named = sorted({f"0{d}" for d in _re.findall(r"\b0([1-4])\b", annotation)})
        box = row.get("positive_box")
        if not named or box not in ("near_label", "far_label"):
            continue                       # nothing to check, or both boxes stained
        vsi = Path(raw_root) / f"Image_{row['tube_id']}.vsi"
        if not vsi.exists():
            continue
        try:
            boxes = assign_boxes(read_vsi_meta(vsi).series)
        except Exception as exc:           # noqa: BLE001 - reported, never fatal
            out.append(f"tube {row['tube_id']}: could not derive boxes ({exc})")
            continue
        target = set(boxes[box])
        if not set(named) <= target:
            out.append(
                f"tube {row['tube_id']}: positive_box={box} resolves to sections "
                f"{sorted(target)} from stage X, but the annotation names {named} "
                f'("{annotation}"). One of the two is wrong -- confirm at the bench.')
    return out


def main(argv):
    csv_path = Path(argv[1]) if len(argv) > 1 else Path(__file__).resolve().parent / "slides.csv"
    print(f"validating {csv_path}")
    if not csv_path.exists():
        print("ERROR: file not found -- this is the one hard blocker for stage 3")
        return 1
    errors, warnings = validate_slides(csv_path)
    for w in warnings:
        print(f"  WARN  {w}")
    for e in errors:
        print(f"  ERROR {e}")
    print(f"{len(errors)} error(s), {len(warnings)} warning(s)")
    if warnings and not errors:
        print("Rows flagged needs_confirmation are NOT usable until confirmed at the bench.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
