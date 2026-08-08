"""Shared helpers for the stage-3 test files (manifest, blinding, QuPath export).

Not a test module -- the name deliberately starts with ``_`` so pytest does not
collect it.  ``tests/test_blinding.py`` contains meta-tests that exercise the
leak detectors defined here, because a *broken detector* is the one failure mode
that would make the blinding tests pass vacuously.

Two problems this file solves:

1.  **Column-name drift.**  The stage-3 modules are being written in parallel to
    a documented function signature, not a documented schema.  Rather than
    hard-coding one spelling of every column, tests ask for a *canonical* name
    and :func:`resolve_column` maps it onto whatever the implementation actually
    emitted, from a short alias list.  A name outside the alias list fails
    loudly with the accepted spellings printed -- it never skips, because a
    silently skipped condition test is exactly as dangerous as a wrong one.

2.  **Naive substring leak-hunting produces false positives.**  Tube IDs are
    two-digit numbers 29-60.  ``"29"`` occurs inside the stage coordinate
    ``12934.5``, ``"32"`` inside the pixel size ``0.325``, ``"38"`` inside the
    width ``16384``.  A blanket ``str(cell).find("29")`` scan would fail on
    perfectly clean data and be switched off within a week.  The scanners below
    separate *exact value* matches (safe on numeric columns) from *token*
    matches inside genuine strings (word-boundary anchored).
"""

from __future__ import annotations

import csv
import math
import random
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SLIDES_CSV = REPO_ROOT / "config" / "slides.csv"

# --------------------------------------------------------------------------
# cohort constants (see CLAUDE_v1.2.md section 4)
# --------------------------------------------------------------------------
ALL_TUBES = tuple(list(range(29, 59)) + [60])

def _repo_root():
    return REPO_ROOT


def _data_root():
    """The raw-data root, from $IHC_DATA_ROOT or config/paths.yaml. None if absent."""
    import os
    env = os.environ.get("IHC_DATA_ROOT")
    if env:
        p = Path(env)
        return p if p.is_dir() else None
    try:
        import yaml
        cfg = yaml.safe_load(open(REPO_ROOT / "config" / "paths.yaml"))
        raw = cfg["roots"]["raw_root"]
        p = Path(str(raw).format(**{k: str(v) for k, v in cfg["roots"].items()
                                    if k != "raw_root"})).expanduser()
        return p if p.is_dir() else None
    except Exception:
        return None



#: Treatment blocks are CONTIGUOUS in tube ID.  This is the whole reason coded
#: IDs must come from a random permutation: any order-preserving scheme
#: reproduces the block structure below exactly.
GROUP_BLOCKS = {
    "Rapamycin Diet": tuple(range(29, 41)),          # 12
    "Extra Control Diet": tuple(range(41, 49)),      # 8
    "Control IP (vehicle)": tuple(range(49, 55)),    # 6
    "Rapamycin IP": (55, 56, 57, 58, 60),            # 5 (59 excluded pre-imaging)
}

#: Tokens that identify a treatment arm.  Any of these appearing in a blinded
#: artefact unblinds it outright.
GROUP_TOKENS = (
    "rapamycin", "vehicle", "extra control", "control diet", "control ip",
    "rapa", "diet",
)

#: Leak vectors named in the spec and verified present in this data.
#: ``stack1`` is the slide LABEL image: printed tube ID plus a DataMatrix
#: barcode.  ``60_`` is tube 60's internal series-name prefix, which QuPath
#: displays in preference to the filename.
FORBIDDEN_SUBSTRINGS = (
    ".vsi", "image_", "stack1", "60_", "_image_", "/", "\\",
)

#: The same list minus the bare path separators, for scanning artefact *files*:
#: a written artefact may legitimately record its own output directory, and
#: failing on that would be noise rather than a finding.
WORK_DIR_FORBIDDEN = (".vsi", "image_", "_image_", "stack1", "60_")

#: Animals whose pixel payload folder is present locally.  The other 23 are
#: index-only, which is a NORMAL state for the manifest, not an error.
def _payload_tubes():
    """Tubes whose pixel payload is on disk right now.

    Derived, not frozen. Payloads arrive in batches over weeks, so a hard-coded list
    turns every legitimate delivery into a red build -- which is how people learn to
    ignore red. Assert *properties* against this (every rescan has an original), never
    membership.
    """
    import glob, os, re
    root = _data_root()
    if root is None:
        return ()
    found = set()
    for p in glob.glob(os.path.join(str(root), "_Image_*_")):
        m = re.search(r"_Image_(\d+)_", os.path.basename(p))
        if m:
            found.add(int(m.group(1)))
    return tuple(sorted(found))


PAYLOAD_TUBES = _payload_tubes()

#: Six animals were scanned with three tissue series, not four.
def _three_section_tubes():
    """Tubes whose .vsi index declares three tissue series rather than four.

    Read from the files. This set genuinely changed when tubes 33, 42 and 54 turned
    out to have a fourth section sitting unscanned in Rescan/.
    """
    import glob, os, re, sys
    root = _data_root()
    if root is None:
        return frozenset()
    sys.path.insert(0, str(_repo_root() / "src"))
    from ihc.ingest.vsi_meta import read_vsi_meta, VsiParseError
    out = set()
    for p in sorted(glob.glob(os.path.join(str(root), "Image_*.vsi"))):
        try:
            meta = read_vsi_meta(p)
        except VsiParseError:
            continue
        if meta.n_tissue_series == 3 and meta.tube_id is not None:
            out.add(meta.tube_id)
    return frozenset(out)


THREE_SECTION_TUBES = _three_section_tubes()

def _needs_confirmation_tubes():
    """Tubes whose ``needs_confirmation`` cell in slides.csv is non-empty.

    Read from the file rather than hard-coded, because this set SHRINKS as the bench
    answers questions, and a frozen constant turns "this has been resolved" into a
    test failure. Tube 37 was the live case until 2026-08-07; it is now empty. The
    mechanism it exercised -- refuse to guess, mark unresolved, exclude from the
    analysis manifest, report loudly -- still needs testing, so the tests read the
    file and adapt to whatever is currently open.
    """
    import csv as _csv
    from pathlib import Path as _Path
    path = _Path(__file__).resolve().parents[1] / "config" / "slides.csv"
    if not path.exists():
        return frozenset()
    with open(path, newline="") as fh:
        return frozenset(
            int(row["tube_id"])
            for row in _csv.DictReader(fh)
            if (row.get("needs_confirmation") or "").strip()
        )


#: Rows in slides.csv whose positive box is contradicted by the bench annotation.
#: Must be excluded from the analysis manifest and reported. Empty as of 2026-08-07.
NEEDS_CONFIRMATION_TUBES = _needs_confirmation_tubes()

#: Slides where BOTH boxes were stained -- these animals have no negative
#: control at all (layouts 4+0 and 3+0).
BOTH_BOXES_TUBES = frozenset({35, 38, 45, 53})


# --------------------------------------------------------------------------
# slides.csv -- the wet-lab record, and the authority on condition
# --------------------------------------------------------------------------
def read_slides_csv(path: Path | str = SLIDES_CSV) -> dict[int, dict[str, str]]:
    """Return ``{tube_id: row}`` from ``config/slides.csv``."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {int(r["tube_id"]): r for r in rows}


def require_cohort_index_files(data_root):
    """Skip unless the whole cohort's ``.vsi`` index files are reachable.

    ``conftest.data_root`` only skips when the directory is *missing*.  A root
    that exists but is empty -- a fresh clone, a wrong ``$IHC_DATA_ROOT``, a
    Dropbox folder that has not synced -- would otherwise let `build_manifest`
    return a near-empty frame and turn every cohort assertion into a failure
    that says nothing about the code.  Partial data is not a weaker test, it is
    a misleading one.
    """
    import pytest

    root = Path(data_root)
    found = sorted(root.glob("Image_*.vsi"))
    if len(found) < len(ALL_TUBES):
        pytest.skip(
            f"only {len(found)} of {len(ALL_TUBES)} Image_*.vsi index files under "
            f"{root}; set $IHC_DATA_ROOT to the full cohort"
        )
    return found


def group_of_tube(tube: int) -> str:
    """Treatment group for ``tube``, from the contiguous blocks above."""
    for name, tubes in GROUP_BLOCKS.items():
        if tube in tubes:
            return name
    raise KeyError(f"tube {tube} is not in the cohort")


# --------------------------------------------------------------------------
# column resolution
# --------------------------------------------------------------------------
#: canonical name -> accepted spellings, most likely first.
COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "tube_id": ("tube_id", "tube", "animal_id", "animal", "mouse_id"),
    "group": ("group", "treatment_group", "treatment"),
    "arm": ("arm", "route", "delivery"),
    "section_label": (
        "section_label", "section", "section_number", "physical_section",
        "section_id",
    ),
    "box": ("box", "box_label", "pap_box", "pap_pen_box", "box_assignment"),
    "condition": (
        "condition", "section_condition", "stain_condition", "staining_condition",
        "positive_negative",
    ),
    "has_negative_control": (
        "has_negative_control", "negative_control_present", "has_negatives",
        "negative_control_available",
    ),
    "payload_present": (
        "payload_present", "has_payload", "pixels_present", "payload_available",
        "payload",
    ),
    "positive_box": ("positive_box", "positive_box_label"),
    "is_rescan": ("is_rescan", "rescan", "rescanned", "is_rescanned"),
    # The categorical form of the same thing: "original" / "rescan".
    "scan_kind": ("scan", "scan_kind", "acquisition", "scan_type"),
    "coded_id": ("coded_id", "code", "blinded_id", "coded_animal_id", "coded_tube"),
    "needs_confirmation": (
        "needs_confirmation", "needs_confirmation_note", "unresolved_reason",
        "confirmation_note",
    ),
    "series_name": ("series_name", "internal_series_name", "name"),
    "series_index": ("series_index", "series", "stack_id", "stack", "series_number"),
    "stage_x_um": ("stage_x_um", "stage_x", "stage_origin_x_um"),
    "acquisition_time": (
        "acquisition_time", "acquisition_timestamp", "acquired_at", "timestamp",
    ),
}


class ColumnNotFound(AssertionError):
    """The manifest has no column matching a canonical name or any alias."""


def find_column(df, canonical: str) -> str | None:
    """Return the actual column name for ``canonical``, or ``None``.

    Matching is case-insensitive and ignores surrounding whitespace, because a
    column called ``Condition`` means the same thing as ``condition`` and a test
    that fails on capitalisation teaches nobody anything.
    """
    aliases = COLUMN_ALIASES.get(canonical, (canonical,))
    lowered = {str(c).strip().lower(): c for c in df.columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    return None


def resolve_column(df, canonical: str) -> str:
    """Like :func:`find_column` but raises with the accepted spellings listed.

    Deliberately raises rather than skipping.  A condition column that cannot be
    located is a schema failure, and skipping would hide it behind a green run.
    """
    found = find_column(df, canonical)
    if found is None:
        raise ColumnNotFound(
            f"the manifest has no {canonical!r} column.\n"
            f"  accepted spellings: {COLUMN_ALIASES.get(canonical, (canonical,))}\n"
            f"  columns present:    {list(df.columns)}"
        )
    return found


def col(df, canonical: str):
    """The resolved column as a Series."""
    return df[resolve_column(df, canonical)]


# --------------------------------------------------------------------------
# value normalisation
# --------------------------------------------------------------------------
_POSITIVE = {"positive", "pos", "p", "+", "primary", "stained", "true"}
_NEGATIVE = {
    "negative", "neg", "n", "-", "control", "negative_control", "no_primary",
    "secondary_only", "false",
}
_UNRESOLVED = {
    "unresolved", "unknown", "pending", "needs_confirmation", "unconfirmed",
    "excluded", "none", "nan", "", "na", "n/a",
}


def norm_condition(value) -> str:
    """Map a condition cell onto ``positive`` / ``negative`` / ``unresolved``.

    An unrecognised value returns ``"?<repr>"`` rather than guessing, so an
    assertion failure names the value that was not understood.

    A null/empty cell counts as ``unresolved``: it is sloppy but it is not
    *wrong*, and the tests that care about positive and negative check for those
    words explicitly, so this leniency cannot make a real failure pass.
    """
    if value is None:
        return "unresolved"
    if isinstance(value, float) and math.isnan(value):
        return "unresolved"
    if isinstance(value, bool):
        return "positive" if value else "negative"
    text = str(value).strip().lower()
    if text in _POSITIVE:
        return "positive"
    if text in _NEGATIVE:
        return "negative"
    if text in _UNRESOLVED:
        return "unresolved"
    return f"?{value!r}"


def norm_box(value) -> str:
    """Map a box cell onto ``near_label`` / ``far_label``, or pass it through."""
    text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if text in {"near_label", "near", "left", "label_end", "low_stage_x"}:
        return "near_label"
    if text in {"far_label", "far", "right", "away_from_label", "high_stage_x"}:
        return "far_label"
    return text


def as_bool(value) -> bool | None:
    """Coerce a truthy/falsy manifest cell to a real bool, or ``None``."""
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, (int,)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "yes", "y", "1", "t"}:
        return True
    if text in {"false", "no", "n", "0", "f"}:
        return False
    return None


def section_key(value) -> str:
    """Normalise a section label to the zero-padded ``"01"`` .. ``"04"`` form."""
    text = str(value).strip()
    if text.endswith(".0"):  # a section number that leaked through as a float
        text = text[:-2]
    digits = re.sub(r"\D", "", text)
    return digits.zfill(2) if digits else text


# --------------------------------------------------------------------------
# row selection
# --------------------------------------------------------------------------
def rescan_mask(df):
    """Boolean Series: True where the row comes from a re-acquisition.

    Two representations are accepted, because both are reasonable and the choice
    is not part of the documented interface: a boolean ``is_rescan`` column, or a
    categorical ``scan`` column reading ``original`` / ``rescan``.  Returns
    ``None`` if the manifest marks rescans in neither way.
    """
    flag = find_column(df, "is_rescan")
    if flag is not None:
        return df[flag].map(as_bool).fillna(False).astype(bool)
    kind = find_column(df, "scan_kind")
    if kind is not None:
        values = df[kind].astype(str).str.strip().str.lower()
        if set(values.unique()) <= {"original", "rescan", "nan", ""}:
            return values.eq("rescan")
    return None


def rows_for_tube(df, tube: int, *, originals_only: bool = True):
    """Every manifest row for ``tube``.

    ``originals_only`` drops rescan rows.  Tubes 51 and 60 were re-acquired at
    standard exposure and those rescans hold two tissue series each (the
    positive box only), so a section-count assertion that forgets to exclude
    them counts 6 rows for a 4-section animal.
    """
    subset = df[col(df, "tube_id").astype("int64") == int(tube)]
    if originals_only:
        mask = rescan_mask(df)
        if mask is not None:
            subset = subset[~mask.loc[subset.index]]
    return subset


def conditions_by_section(df, tube: int) -> dict[str, str]:
    """``{"01": "positive", ...}`` for one animal, from the original scan."""
    subset = rows_for_tube(df, tube)
    labels = col(subset, "section_label").map(section_key)
    conditions = col(subset, "condition").map(norm_condition)
    return dict(zip(labels, conditions))


def boxes_by_section(df, tube: int) -> dict[str, str]:
    """``{"01": "far_label", ...}`` for one animal, from the original scan."""
    subset = rows_for_tube(df, tube)
    labels = col(subset, "section_label").map(section_key)
    boxes = col(subset, "box").map(norm_box)
    return dict(zip(labels, boxes))


# --------------------------------------------------------------------------
# leak scanning
# --------------------------------------------------------------------------
def iter_cells(df):
    """Yield ``(column_name, row_position, value)`` for every cell."""
    for column in df.columns:
        for position, value in enumerate(df[column].tolist()):
            yield str(column), position, value


def _tube_value_forms(tube: int) -> set:
    """Every literal form a tube ID could take as a whole cell value."""
    return {tube, float(tube), str(tube), f"{tube:02d}", f"{tube}.0", f"Image_{tube}"}


def exact_tube_id_hits(df, tubes=ALL_TUBES) -> list[tuple[str, int, object, int]]:
    """Cells whose *whole value* is a tube ID.

    Safe to run over numeric columns: no legitimate measurement in this dataset
    (stage coordinates in um, pixel size 0.325, widths ~16000, exposures in ms)
    is exactly equal to an integer in 29..60.  A section count of 3 or 4 and a
    series index are both outside the range too.
    """
    forms = {tube: _tube_value_forms(tube) for tube in tubes}
    hits = []
    for column, position, value in iter_cells(df):
        if isinstance(value, str):
            probe = value.strip()
        else:
            probe = value
        for tube, candidates in forms.items():
            if probe in candidates:
                hits.append((column, position, value, tube))
    return hits


def _token_pattern(tubes=ALL_TUBES) -> re.Pattern:
    """Pattern matching a tube ID as a standalone number.

    The four lookarounds are all load-bearing, and each was added because the
    naive version produced a wrong answer on this cohort's real values:

    * ``(?<!\\d)`` / ``(?!\\d)`` -- so ``16384`` is not read as tube 38 and
      ``1840.0`` is not read as tube 40.
    * ``(?<!\\d\\.)`` / ``(?!\\.\\d)`` -- so the *decimal* ``0.40`` and ``40.0``
      are not read as tube 40 either.  Note this still matches ``Image_29.vsi``,
      because the character after the dot is a letter, not a digit.
    """
    alternatives = "|".join(str(t) for t in sorted(tubes, reverse=True))
    return re.compile(rf"(?<!\d)(?<!\d\.)(?:{alternatives})(?!\d)(?!\.\d)")


_TUBE_TOKEN_RE = _token_pattern()


def text_tube_id_hits(text: str, tubes=ALL_TUBES) -> list[str]:
    """Tube IDs appearing as standalone numbers anywhere in ``text``.

    Used to scan artefact *files* -- a blinded CSV that never mentions a tube ID
    in a cell can still leak one in a header, a comment or a provenance line.
    """
    wanted = {str(t) for t in tubes}
    return sorted({m.group(0) for m in _TUBE_TOKEN_RE.finditer(text) if m.group(0) in wanted})


#: Filename stems that carry the tube ID in this cohort.
_FILENAME_TUBE_RE = re.compile(r"_?image[_ ]?(\d{2})_?", re.IGNORECASE)

#: A number is an animal identifier when an identifying word IMMEDIATELY
#: precedes it -- "tube 49", "animal_49", "Slide: 49".  Mere co-occurrence on
#: the line is not enough: `"n_animals": 31` is the cohort size and "29 of 31
#: slides" is a fact from the published spec, and neither names a mouse.
_IDENTIFYING_LABEL_RE = re.compile(
    r"(?:tube|animal|mouse|slide|specimen)s?[\s_:=#.-]*(\d{2})(?!\d)",
    re.IGNORECASE,
)


def associated_tube_id_hits(text: str, codes=None, tubes=ALL_TUBES) -> list[str]:
    """Tube IDs in ``text`` that are *attached to something*, line by line.

    Why not just scan for the digits: a legitimate artefact reports counts.
    ``"n_animals": 31`` is the cohort size, ``"29 of 31 slides"`` is a fact from
    the published spec, and neither tells anyone which mouse is which.  A blanket
    digit scan flags them, the test gets muted, and the real leak walks through.

    What actually constitutes a leak is an *association*: a tube ID next to a
    code (that is the key), a tube ID next to the word "tube"/"animal"/"slide"
    (that is a label), or a filename stem such as ``Image_49`` / ``_Image_49_``
    (that is a path).  Those are what this reports.

    Args:
        text: the file contents.
        codes: ``{tube: code}``.  When given, any line holding both a code and a
            tube ID is reported -- that single line reconstructs the key.
    """
    hits: set[str] = []  # type: ignore[assignment]
    hits = set()
    wanted = {str(t) for t in tubes}
    code_values = {str(c) for c in (codes or {}).values()}

    for stem_match in _FILENAME_TUBE_RE.finditer(text):
        if stem_match.group(1) in wanted:
            hits.add(f"filename stem {stem_match.group(0)!r}")

    for line in text.splitlines():
        found = {m.group(0) for m in _TUBE_TOKEN_RE.finditer(line)} & wanted
        if not found:
            continue
        if code_values and any(code in line for code in code_values):
            hits.add(f"tube {sorted(found)} on the same line as a code: {line.strip()[:120]!r}")
            continue
        labelled = {m.group(1) for m in _IDENTIFYING_LABEL_RE.finditer(line)} & found
        if labelled:
            hits.add(f"tube {sorted(labelled)} labelled as an animal: {line.strip()[:120]!r}")
    return sorted(hits)


def text_forbidden_hits(text: str, tokens=None) -> list[str]:
    """Forbidden tokens (group names, path fragments, leak markers) in ``text``."""
    if tokens is None:
        tokens = FORBIDDEN_SUBSTRINGS + GROUP_TOKENS
    lowered = text.lower()
    return sorted({t for t in tokens if t in lowered})


def substring_tube_id_hits(df, tubes=ALL_TUBES) -> list[tuple[str, int, str, str]]:
    """Tube IDs appearing *inside string cells*, anchored to digit boundaries.

    Only genuine ``str`` cells are scanned.  Scanning stringified floats would
    flag ``0.325`` for tube 32 and ``16384`` for tube 38, which is noise.  The
    digit-boundary anchor means a random code such as ``"SLK2291"`` is not
    flagged, while ``"slide-29"``, ``"Image_29.vsi"`` and an ISO date ending
    ``-07-29`` are.  An acquisition timestamp being flagged is a true positive,
    not a nuisance: the cohort was scanned in ascending tube order, so the
    timestamp is a perfect group proxy.
    """
    wanted = {str(t) for t in tubes}
    hits = []
    for column, position, value in iter_cells(df):
        if not isinstance(value, str):
            continue
        for match in _TUBE_TOKEN_RE.finditer(value):
            if match.group(0) in wanted:
                hits.append((column, position, value, match.group(0)))
    return hits


def forbidden_substring_hits(
    df, tokens=FORBIDDEN_SUBSTRINGS + GROUP_TOKENS
) -> list[tuple[str, int, str, str]]:
    """String cells containing a group name, a path fragment or a known leak token."""
    hits = []
    for column, position, value in iter_cells(df):
        if not isinstance(value, str):
            continue
        lowered = value.lower()
        for token in tokens:
            if token in lowered:
                hits.append((column, position, value, token))
    return hits


def column_name_hits(df, tokens=FORBIDDEN_SUBSTRINGS + GROUP_TOKENS) -> list[str]:
    """Column *names* that themselves carry a forbidden token."""
    out = []
    for column in df.columns:
        lowered = str(column).lower()
        for token in tokens:
            if token in lowered and token not in {"/", "\\"}:
                out.append(f"{column} (contains {token!r})")
                break
    return out


# --------------------------------------------------------------------------
# partitions -- "does this column separate the treatment groups?"
# --------------------------------------------------------------------------
def partition_from(values) -> frozenset:
    """Partition of row positions induced by a sequence of values.

    ``None`` and ``NaN`` collapse to a single ``"<null>"`` bucket so that two
    missing values are treated as equal, which is what a human reading the
    column would assume.
    """
    buckets: dict[str, set[int]] = {}
    for position, value in enumerate(values):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            key = "<null>"
        else:
            key = repr(value)
        buckets.setdefault(key, set()).add(position)
    return frozenset(frozenset(members) for members in buckets.values())


def partition_of_groups(groups) -> frozenset:
    """The true treatment-group partition of the same row positions."""
    return partition_from(list(groups))


# --------------------------------------------------------------------------
# block structure -- the statistic that catches order-preserving codes
# --------------------------------------------------------------------------
def same_group_adjacent_pairs(group_sequence) -> int:
    """Adjacent pairs sharing a treatment group, in the given order.

    Tube IDs run in contiguous blocks of 12, 8, 6 and 5, so tube order gives 27
    of a possible 30 same-group adjacencies.  A random permutation gives about
    7.7 on average.  The gap is what makes this statistic decisive.
    """
    sequence = list(group_sequence)
    return sum(1 for a, b in zip(sequence, sequence[1:]) if a == b)


def adjacency_null(group_sequence, *, draws: int = 20000, seed: int = 20260807):
    """Monte-Carlo null distribution of :func:`same_group_adjacent_pairs`.

    Returns a sorted list.  The null is generated with a *fixed* seed so the
    whole test is deterministic: given a fixed blinding seed, the test either
    always passes or always fails, never intermittently.
    """
    labels = list(group_sequence)
    rng = random.Random(seed)
    out = []
    for _ in range(draws):
        rng.shuffle(labels)
        out.append(same_group_adjacent_pairs(labels))
    out.sort()
    return out


def quantile(sorted_values, q: float):
    """The ``q``-quantile of an already-sorted list (nearest-rank)."""
    if not sorted_values:
        raise ValueError("empty distribution")
    rank = max(1, math.ceil(q * len(sorted_values)))
    return sorted_values[min(rank, len(sorted_values)) - 1]


def mean(values) -> float:
    values = list(values)
    return sum(values) / len(values)


def stdev(values) -> float:
    values = list(values)
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((v - mu) ** 2 for v in values) / (len(values) - 1))


# --------------------------------------------------------------------------
# blinded-row -> tube mapping
# --------------------------------------------------------------------------
def code_column(blinded_df, codes: dict) -> str:
    """The blinded column holding coded IDs.

    Found by *value*, not by name: the column whose non-null values are all
    members of ``codes.values()``.  This keeps the leak tests independent of
    what the implementation decided to call it.
    """
    wanted = {str(v) for v in codes.values()}
    for column in blinded_df.columns:
        values = [v for v in blinded_df[column].tolist() if v is not None]
        values = [v for v in values if not (isinstance(v, float) and math.isnan(v))]
        if not values:
            continue
        if all(str(v) in wanted for v in values):
            return str(column)
    raise AssertionError(
        "no column in the blinded manifest holds coded IDs.\n"
        f"  codes issued: {sorted(wanted)[:5]}...\n"
        f"  columns:      {list(blinded_df.columns)}"
    )


def tubes_of_blinded_rows(blinded_df, codes: dict) -> list[int]:
    """The true tube ID behind each blinded row, via the code column."""
    inverse = {str(code): tube for tube, code in codes.items()}
    column = code_column(blinded_df, codes)
    return [inverse[str(value)] for value in blinded_df[column].tolist()]


# --------------------------------------------------------------------------
# synthetic manifests (used where real data would be slower and no clearer)
# --------------------------------------------------------------------------
def synthetic_private_manifest(pd, tubes=ALL_TUBES, *, condition_of=None):
    """A private manifest with the columns the blinding stage consumes.

    Deliberately includes every known leak vector -- tube ID in a path, the tube
    ID inside the internal series name (``60_20x_...``, which QuPath displays in
    preference to the filename), a monotone acquisition timestamp, and the two
    non-standard exposures -- so that ``split_manifest`` has something real to
    strip.  A blinding test run against a manifest with nothing to leak proves
    nothing.
    """
    slides = read_slides_csv()
    rows = []
    clock = 0
    for tube in tubes:  # ascending: the real scanning order, and the leak
        record = slides.get(tube, {})
        n_sections = 3 if tube in THREE_SECTION_TUBES else 4
        positive_box = record.get("positive_box", "near_label")
        for index in range(n_sections):
            label = f"{index + 1:02d}"
            box = "near_label" if index < 2 else "far_label"
            if positive_box == "both":
                condition = "positive"
            else:
                condition = "positive" if box == positive_box else "negative"
            if condition_of is not None:
                condition = condition_of(tube, label, condition)
            clock += 137
            rows.append(
                {
                    "tube_id": tube,
                    "group": group_of_tube(tube),
                    "arm": "diet" if tube <= 48 else "ip",
                    "section_label": label,
                    "box": box,
                    "condition": condition,
                    "positive_box": positive_box,
                    "has_negative_control": positive_box != "both",
                    "payload_present": tube in PAYLOAD_TUBES,
                    "is_rescan": False,
                    "series_name": (
                        f"{tube}_20x_DAPI, FITC, Cy3_{label}"
                        if tube == 60
                        else f"20x_DAPI, FITC, Cy3_{label}"
                    ),
                    "series_index": 10002 + 3 * index,
                    "original_path": f"/RawData/Image_{tube}.vsi",
                    "stage_x_um": 12000.0 + 4200.0 * index,
                    "pixel_size_um": 0.325,
                    "width_px": 16384,
                    "height_px": 20480,
                    "exposure_dapi_ms": 60.52 if tube == 51 else 128.55,
                    "exposure_fitc_ms": 240.82 if tube == 51 else 397.93,
                    "exposure_cy3_ms": (
                        145.68 if tube == 51 else 397.93 if tube == 60 else 1840.0
                    ),
                    "acquisition_time": f"2026-07-29T{10 + clock // 3600:02d}:"
                    f"{(clock // 60) % 60:02d}:{clock % 60:02d}+00:00",
                    "acquisition_order": clock,
                }
            )
    return pd.DataFrame(rows)
