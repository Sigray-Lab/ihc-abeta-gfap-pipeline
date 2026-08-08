#!/usr/bin/env python3
"""Split the ingest manifest into a private provenance record and a blinded one.

Why this module exists
----------------------
The person who draws the anatomical regions must not know which animal is which,
and in this cohort that is harder than it sounds. wet-lab performed the staining and
the imaging, so he can infer treatment group from a tube ID alone; and the tube ID
is *everywhere*. Every one of the following was verified present in this data:

* the file path and the file name contain the tube ID (``Image_29.vsi``);
* ``stack1`` is the slide LABEL image, which shows ``1007344 - 29`` as printed text
  **and** as a DataMatrix barcode encoding the same string;
* the internal series names carry it too -- tube 60's are
  ``60_20x_DAPI, FITC, Cy3_01`` -- and QuPath displays *internal* names, not file
  names, so renaming the file achieves nothing;
* VSI tags 2061 and 120635 hold it as a string;
* the acquisition timestamp is a perfect proxy: scanning ran in ascending tube
  order, and tube order is ascending group order;
* exposure identifies two animals outright: tubes 51 and 60 are the only slides
  acquired at non-standard exposure.

And tube IDs run in contiguous treatment blocks --  29-40 Rapamycin Diet, 41-48
Extra Control Diet, 49-54 Control IP, 55-60 Rapamycin IP -- so **any order-preserving
coding scheme reproduces the group structure exactly**. Sequential codes in file
order, a hash, ``tube_id * 7 + 3``: all of them keep the blocks contiguous and all of
them are worthless. Codes therefore come from a random permutation with a recorded
seed (spec §2, ADR-0009), and :func:`audit_blinded` tests that they did.

What the four public functions do
---------------------------------
:func:`generate_codes`
    ``tube_id -> code`` from ``random.Random(seed).shuffle`` over a fixed label pool.
:func:`split_manifest`
    one manifest in, two out: everything for the custodian, a safe subset for the
    analyst. Refuses outright if any row's condition is unresolved.
:func:`write_blinded`
    writes the key + the private manifest to a mode-0700 custodian directory and the
    blinded manifest to the work directory, and refuses to write the blinded manifest
    at all if the audit finds a hard leak.
:func:`audit_blinded`
    tries to break the blinding and returns everything it finds.

The severity ladder used by :func:`audit_blinded`
-------------------------------------------------
Every finding is a string starting with one of two tokens. An empty list means
nothing was found.

``LEAK:``
    The blinded artefact **directly** exposes identity: a tube ID, a group or arm
    label, a path, a column that is one-to-one with animal (and therefore joinable
    against any un-blinded table), a column that is monotone in tube ID, or a code
    order that preserves the treatment blocks. These are hard stops -- ``write_blinded``
    refuses to write the blinded manifest.

``RISK:``
    The blinded artefact **narrows** identity without exposing it: a value carried by
    only one or two animals, which anyone holding the bench record can resolve to a
    named animal. Reported, recorded in the report JSON, printed by ``./ihc blind`` --
    but it needs a human judgement rather than a machine veto, because some of these
    are unavoidable properties of the cohort.

That distinction is load-bearing here, because the obvious "safe" replacement for raw
exposure -- a boolean ``exposure_is_standard`` -- is itself a ``RISK`` finding: exactly
two of thirty-one animals are ``False``, the spec, the execution plan and ADR-0006 all
name those two animals in writing, and one of them is Control IP while the other is
Rapamycin IP. Coarsening the value does not remove the leak, it only makes it two bits
instead of a float. Use the rescans (acquired at standard exposure) if you want it gone.

Public API
----------
:func:`generate_codes`, :func:`split_manifest`, :func:`write_blinded`,
:func:`audit_blinded`; exceptions :class:`UnresolvedConditionError`,
:class:`BlindingLeakError`, :class:`CustodianPathError`.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import os
import random
import re
import statistics
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence

__all__ = [
    "generate_codes",
    "normalise_seed",
    "split_manifest",
    "write_blinded",
    "audit_blinded",
    "UnresolvedConditionError",
    "BlindingLeakError",
    "CustodianPathError",
    "CODE_POOL_LETTERS",
    "CODE_POOL_NUMBERS",
    "BLINDED_COLUMNS",
    "NEVER_BLINDED",
]


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

#: Letters used to build code labels. ``I`` and ``O`` are omitted because they are
#: read as ``1`` and ``0`` on a printout, and a code that is transcribed wrongly by
#: hand cannot be un-blinded later.
CODE_POOL_LETTERS = "ABCDEFGHJKLMNPQRSTUVWXYZ"

#: Numeric part of a code label. Kept strictly below the lowest tube ID in this
#: cohort (29) so that a code can never coincidentally read as a tube number.
CODE_POOL_NUMBERS = range(1, 26)

#: The pool must be this many times larger than the cohort. A code drawn from a pool
#: several times its own size tells an observer nothing about how many animals there
#: are, nor where any animal sits in the sequence.
CODE_POOL_HEADROOM = 4

#: Columns allowed through to the blinded manifest, by exact name. Anything not named
#: here is **dropped and reported** -- the default is to withhold, because a column
#: nobody thought about is exactly the one that leaks. Widen it deliberately via
#: ``split_manifest(..., extra_blinded_columns=[...])``, and read the audit afterwards.
BLINDED_COLUMNS = (
    "code",                   # the coded ID; the only identifier that crosses over
    "row_kind",               # "section" -- schema marker, constant
    "analysis_include",       # constant True in an analysis manifest; kept for schema
    # The PHYSICAL section identity, derived from stage X -- not the scanner's label.
    # These differ, and using the scanner's label here was an actual bug: tube 51's
    # rescan numbers its two sections 01/02, but physically they are the original
    # slide's 03/04. Keying on the scanner label collapsed the original's 01
    # (negative control) and the rescan's 01 (positive) onto the same (code, label)
    # with OPPOSITE conditions, so a reader of the blinded manifest could not tell
    # which row was the control. Physical identity is the only stable key across scans.
    "physical_section_label",
    "section_label",          # kept for reference; scanner acquisition order within a scan
    "replicate",              # 1..n within a section; see split_manifest
    "series_index",           # which series inside the file, for opening the right image
    "stack_id",               # 10002/10005/10008/10011 -- same job, manifest.py's name
    "box",                    # near_label / far_label, from stage X
    "condition",              # positive / negative -- from slides.csv, never from pixels
    "region",                 # hippocampus / isocortex, if the manifest carries it
    "marker",                 # abeta / gfap, if the manifest carries it
    "channel",                # DAPI / FITC / Cy3
    "n_sections_on_slide",    # 3 or 4
    "n_sections_recorded",
    "n_positive_sections_on_slide",
    "has_negative_control",   # false for the four 4+0 slides
    "pixel_size_um",          # rounded; see PIXEL_SIZE_DECIMALS
    "exposure_is_standard",   # derived boolean; NOT the raw exposure. See the docstring.
)

#: Columns that must never reach the blinded manifest whatever the allow-list says,
#: matched as case-insensitive substrings of the column name. Each one is a verified
#: leak vector in this dataset, not a hypothetical.
NEVER_BLINDED = (
    "tube",            # the ID itself, under any name
    "animal",          # ... and its synonyms
    "subject",
    "mouse",
    "group",           # treatment group
    "arm",             # diet / ip
    "treatment",
    "path",            # any path is a file name is a tube ID
    "file",
    "dir",
    "folder",
    "uri",
    "url",
    "series_name",     # tube 60's internal series names begin "60_"
    "acquisition",     # scan order == tube order == group order
    "acq",             # ... and its abbreviation
    "scan",            # scan_day, scan_order, rescan
    "timestamp",
    "time",            # a blinded manifest has no legitimate time column
    "date",
    "datetime",
    "batch",
    "session",
    "exposure_ms",     # raw exposure identifies tubes 51 and 60
    "exposure_dapi",
    "exposure_fitc",
    "exposure_cy3",
    "exposure_us",
    "sha",             # a checksum is a per-file fingerprint
    "checksum",
    "hash",
    "md5",
    "label_image",     # stack1 shows the tube ID as text and as a DataMatrix barcode
    "seed",            # the key must never travel with the blinded data
    "code_seed",
    "depth_index",     # slide depth index isolates single animals (values 4 and 6 are unique)
    "rescan",          # only tubes 51 and 60 (and four repeat-scanned slides) are rescans
    "needs_confirmation",
    "annotation",      # free text from the bench record; names animals and sections
    "warning",         # warnings quote series names and paths
)

#: Values of ``condition`` that mean "we do not know", any one of which stops the run.
#: Tube 37 is the live case: ``slides.csv`` names ``far_label`` while the bench
#: annotation names section 01, which is in ``near_label``. Blinding that row would
#: bury a contradiction under a code, where nobody would ever find it again.
UNRESOLVED_CONDITIONS = frozenset(
    {"", "unresolved", "unknown", "pending", "PENDING_PI_DECISION", "nan",
     "none", "tbd", "needs_confirmation", "?"}
)

#: Pixel size is 0.32500-0.32502 um per series -- enough precision to fingerprint an
#: individual series. Rounded to this many decimals it is a constant, which is what
#: downstream actually needs.
PIXEL_SIZE_DECIMALS = 4

#: Exposure triples are compared to the cohort mode within this relative tolerance.
EXPOSURE_RTOL = 1e-3

#: Monte Carlo settings for the code-order tests. The seed is a fixed constant and is
#: deliberately NOT the blinding seed: the audit must be reproducible by someone who
#: does not hold the key, and must not vary with it.
AUDIT_MC_PERMUTATIONS = 20000
AUDIT_MC_SEED = 1000003          # a prime, deliberately not date-shaped: it is not a key

#: |rho| at or above this between tube ID and code rank is a hard fail on its own,
#: independent of the Monte Carlo p-value.
MAX_ABS_SPEARMAN = 0.5

#: A column whose per-animal values are monotone in tube ID at this |rho| or above is
#: an order-preserving identifier -- the acquisition-timestamp case.
COLUMN_MONOTONE_RHO = 0.9

#: p-value below which a Monte Carlo tail counts as a failure.
MC_ALPHA = 0.01

#: p-value below which a code order is not a failure but is worth re-drawing. Choosing
#: another seed when the draw looks structured is restricted randomisation, not cheating.
MC_ADVISORY_ALPHA = 0.05

#: A value carried by this many animals or fewer singles them out.
RARE_VALUE_MAX_ANIMALS = 2

_PENDING = "PENDING_PI_DECISION"

_KEY_CSV_NAME = "blinding_key.csv"
_KEY_JSON_NAME = "blinding_key.json"
_PROVENANCE_CSV_NAME = "provenance_manifest.csv"
_BLINDED_CSV_NAME = "analysis_manifest_blinded.csv"
_REPORT_JSON_NAME = "blinding_report.json"

#: Generic treatment vocabulary, checked in addition to whatever labels the private
#: manifest actually contains, so that a paraphrase does not slip through.
_TREATMENT_WORDS = (
    "rapamycin", "sirolimus", "vehicle", "saline", "placebo", "untreated",
    "diet", "chow", "intraperitoneal",
)

#: Phrases that legitimately contain the word "control" in a blinded manifest. They
#: are masked out before the group-label search, so ``has_negative_control`` does not
#: read as the group label "Control IP (vehicle)".
_BENIGN_CONTROL = re.compile(
    r"(?:\bno[_\- ]?|\bhas[_\- ]?)?(?:negative|positive|neg|pos)[_\- ]?controls?\b",
    re.IGNORECASE,
)

#: Numbers that are structurally *not* tube IDs, masked out before the tube-ID search
#: so that the minute field of ``2026-07-18T23:30:42+00:00`` does not read as tube 30
#: and an exposure of ``145.677`` does not read as tube 45. Timestamps and decimals are
#: still caught, by the uniqueness, monotonicity and forbidden-name tests -- and those
#: catch them for the right reason, which matters when a human has to act on the report.
_NUMERIC_NOISE = re.compile(
    r"\d{4}-\d{2}-\d{2}(?:[T ]\d{2}:\d{2}:\d{2}(?:[.,]\d+)?(?:Z|[+-]\d{2}:?\d{2})?)?"
    r"|\d{1,2}:\d{2}(?::\d{2})?"
    r"|\d+[.,]\d+"
)

_PATHISH = re.compile(
    r"(?:^|[\s,;=])(?:/|~/|\.\./|\./)"          # posix absolute or relative path
    r"|[A-Za-z]:[\\/]"                           # windows drive
    r"|\\\\[A-Za-z0-9._-]+\\"                    # UNC share
    r"|\.(?:vsi|ets|tif|tiff|ome|czi|ndpi|svs|png|jpe?g|json|csv|xlsx?|qpproj|qpdata)\b"
    r"|\bImage_\d+\b"                            # the raw file stem in this cohort
    r"|\b_Image_\d+_\b",                         # the payload folder name
    re.IGNORECASE,
)


# --------------------------------------------------------------------------- #
# Exceptions
# --------------------------------------------------------------------------- #


class UnresolvedConditionError(ValueError):
    """A row's positive/negative condition is not resolved, so nothing may be blinded.

    Raised by :func:`split_manifest`. Condition is read from ``config/slides.csv`` and
    never from pixels (spec §2, ADR-0003). A row whose condition is unresolved, or
    whose ``needs_confirmation`` cell is non-empty, must be excluded from the analysis
    manifest and *reported* -- not coded. Once it has a code, the contradiction is
    invisible: the numbers it produces will look entirely plausible.
    """


class BlindingLeakError(RuntimeError):
    """The audit found a hard leak, so the blinded manifest was not written."""


class CustodianPathError(ValueError):
    """The custodian directory is somewhere the key must never live."""


# --------------------------------------------------------------------------- #
# Codes
# --------------------------------------------------------------------------- #


def normalise_seed(seed):
    """Canonical form of a seed, so the CLI and the API cannot disagree.

    ``random.Random(42)`` and ``random.Random("42")`` are *different* generators and
    produce *different* permutations. ``argparse`` hands over a string, a test hands
    over an int, and the two would silently produce two different keys for the same
    recorded seed -- which is not a hypothetical, it is the obvious way to lose a
    cohort's coding. A seed that reads as an integer is therefore always used as an
    integer; anything else is used as its string form.
    """
    if isinstance(seed, bool):                      # bool is an int; refuse the ambiguity
        raise ValueError("a boolean is not a seed")
    if isinstance(seed, int):
        return seed
    text = str(seed).strip()
    try:
        return int(text)
    except ValueError:
        return text


def _code_pool() -> list[str]:
    """Every candidate code label, in a fixed order that does not depend on the seed."""
    return [f"{letter}{number:02d}"
            for letter in CODE_POOL_LETTERS
            for number in CODE_POOL_NUMBERS]


#: Minimum digits for a numeric seed. `secrets.randbits(128)` gives ~39.
MIN_NUMERIC_SEED_DIGITS = 12

_DATE_SHAPED = re.compile(r"^(19|20)\d{2}(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])$")


def _reject_guessable_seed(seed) -> None:
    """Refuse a seed small enough to brute-force.

    The mapping is deterministic in ``(sorted tube list, seed)``. The cohort roster is
    public -- it is a constant in this repository -- and the issued codes are in the
    blinded manifest that the analyst holds. So an attacker does not need to break
    anything: they enumerate seeds, generate the mapping, and compare against the codes
    they already have. A date is ~365 tries. A two-digit number is 100.

    This exists because an earlier version of this project claimed such a guard in its
    decision log while the code had none -- a false assurance, which is worse than a
    known gap, because nobody re-checks a documented control.

    Draw seeds with ``secrets.randbits(128)``.
    """
    text = str(seed).strip()
    if _DATE_SHAPED.match(text):
        raise ValueError(
            f"seed {text!r} is date-shaped (YYYYMMDD). The code mapping is reproducible "
            f"from (cohort, seed), the cohort roster is public in this repository, and "
            f"the codes are in the blinded manifest -- so a date is about 365 guesses. "
            f"Use secrets.randbits(128).")
    if text.isdigit():
        if len(text) < MIN_NUMERIC_SEED_DIGITS:
            raise ValueError(
                f"seed {text!r} has {len(text)} digits; at least "
                f"{MIN_NUMERIC_SEED_DIGITS} are required. A short numeric seed is "
                f"enumerable in seconds. Use secrets.randbits(128).")
        return
    # Non-numeric seeds are allowed but must not be a guessable word or phrase. There is
    # no way to test "unguessable", so require enough length and character variety that
    # a wordlist attack is not trivial.
    if len(text) < 16 or len(set(text)) < 8:
        raise ValueError(
            f"seed {text!r} is too simple: a non-numeric seed needs at least 16 "
            f"characters and 8 distinct ones, or it is a wordlist entry. Use "
            f"secrets.randbits(128) instead of inventing one.")


def generate_codes(tube_ids: Iterable[int], *, seed, existing: dict | None = None,
                   allow_weak_seed: bool = False) -> dict:
    """Map each tube ID to a code label drawn from a random permutation.

    The permutation is the whole point. Tube IDs run in contiguous treatment blocks,
    so a code that preserves order -- sequential in file-iteration order, a hash, an
    arithmetic transform -- reproduces the group structure exactly and blinds nobody.
    ``random.Random(seed).shuffle`` over a fixed label pool, zipped to the *sorted*
    tube IDs, makes the mapping reproducible from the seed alone and unpredictable
    without it.

    Labels look like ``A07`` or ``T19``: they carry no ordering anybody would read
    into them, they are short enough to appear in a QuPath project name, and their
    numeric part stops at 25 so a code can never be mistaken for a tube number
    (the lowest tube in this cohort is 29). The pool is
    ``CODE_POOL_HEADROOM`` times larger than the cohort, so the set of issued codes
    does not reveal the cohort size either.

    **The mapping is stable under a re-run, but NOT under a change of cohort.** Same
    seed and same animals gives the identical mapping, and the input order does not
    matter. But the codes are assigned by zipping a shuffled pool to the sorted tube
    list, so adding or removing one animal shifts almost every other assignment --
    measured: dropping a single tube left only 8 of 30 animals with their original
    code. That would orphan every measurement already made.

    ``existing`` closes that hole structurally rather than by asking people to be
    careful. Pass the current key and codes already issued are kept verbatim; only
    animals not yet in it draw from the remaining pool. The key becomes append-only,
    so a late-arriving animal can never rescramble the cohort.

    Args:
        tube_ids: The animals to code. Duplicates are an error.
        seed: Anything ``random.Random`` accepts -- an int is what the custodian
            should record. Keyword-only, with no default, deliberately: a default
            seed is a published key.
        existing: An already-issued ``{tube_id: code}``. Those assignments are
            preserved exactly; new animals get codes not already in use.

    Returns:
        ``{tube_id: code}``, ordered by tube ID.

    Raises:
        ValueError: empty input, duplicate tube IDs, no seed, the
            ``PENDING_PI_DECISION`` sentinel as a seed, or a pool too small for the
            cohort.

    Example:
        >>> generate_codes([29, 30, 31], seed=1)["29"] if False else None
    """
    tubes = [int(t) for t in tube_ids]
    if not tubes:
        raise ValueError("no tube IDs to code")
    if len(set(tubes)) != len(tubes):
        dupes = sorted({t for t in tubes if tubes.count(t) > 1})
        raise ValueError(f"duplicate tube IDs: {dupes}")
    if seed is None or seed == "":
        raise ValueError(
            "generate_codes() needs an explicit seed. It is the key: without it the "
            "mapping cannot be reproduced, and with a default it is not a key at all.")
    if seed == _PENDING:
        raise ValueError(
            f"the blinding seed is still {_PENDING!r} in config/config.yaml "
            "(decision D-14). The custodian chooses it and records it outside git; "
            "it is never committed.")

    # allow_weak_seed exists for tests, which legitimately need small deterministic seeds
    # to check that the same seed reproduces a mapping and a different one does not. It is
    # keyword-only, defaults to False, and no production call site passes it -- the CLI
    # never does. Anything that reaches a real key goes through the guard.
    if not allow_weak_seed:
        _reject_guessable_seed(seed)

    tubes.sort()
    seed = normalise_seed(seed)
    pool = _code_pool()
    if len(pool) < CODE_POOL_HEADROOM * len(tubes):
        raise ValueError(
            f"code pool holds {len(pool)} labels, which is less than "
            f"{CODE_POOL_HEADROOM}x the {len(tubes)} animals. Widen "
            f"CODE_POOL_LETTERS or CODE_POOL_NUMBERS.")

    random.Random(seed).shuffle(pool)

    if existing:
        # Append-only. Keep every code already issued and draw for the rest from what
        # is left. Without this, adding one late-arriving animal reshuffles almost the
        # whole cohort -- measured: dropping a single tube left only 8 of 30 codes
        # unchanged -- and every measurement made under the old codes is orphaned.
        existing = {int(t): str(c) for t, c in existing.items()}
        taken = set(existing.values())
        kept = {t: existing[t] for t in tubes if t in existing}
        remaining = [label for label in pool if label not in taken]
        need = [t for t in tubes if t not in existing]
        if len(remaining) < len(need):
            raise ValueError(
                f"{len(need)} new animal(s) to code but only {len(remaining)} unused "
                f"label(s) left in the pool. Widen CODE_POOL_LETTERS or "
                f"CODE_POOL_NUMBERS.")
        kept.update(dict(zip(need, remaining)))
        codes = {t: kept[t] for t in tubes}
    else:
        codes = dict(zip(tubes, pool[: len(tubes)]))

    # A code that contains its own animal's tube number would undo the whole exercise.
    # Impossible with the current pool (numbers stop at 25, tubes start at 29), so this
    # is a guard against a later widening of the pool rather than a live risk.
    for tube, code in codes.items():
        if re.search(rf"(?<!\d){tube}(?!\d)", code):
            raise ValueError(
                f"code {code!r} contains tube ID {tube}; change the seed or the pool")
    return codes


# --------------------------------------------------------------------------- #
# Splitting
# --------------------------------------------------------------------------- #


def _pandas():
    import pandas as pd                                       # noqa: PLC0415
    return pd


def _is_missing(value) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def _text(value) -> str:
    return "" if _is_missing(value) else str(value)


def _column(df, *names):
    """First column present out of `names`, case-insensitively. None if absent."""
    lookup = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    return None


def _unresolved_rows(df):
    """Row indices whose condition is unresolved or which are flagged for confirmation."""
    bad = []
    condition_col = _column(df, "condition")
    confirm_col = _column(df, "needs_confirmation")
    for idx, row in df.iterrows():
        reasons = []
        if condition_col is not None:
            value = _text(row[condition_col]).strip()
            if value.lower() in UNRESOLVED_CONDITIONS or value == _PENDING:
                reasons.append(f"condition={value!r}")
        if confirm_col is not None:
            note = _text(row[confirm_col]).strip()
            # An empty needs_confirmation cell survives a CSV round-trip as "", as NaN,
            # or -- if something inferred the empty column's dtype -- as False or 0.
            # None of those is a flag, and treating them as one would refuse the whole
            # cohort, which is a failure mode that looks exactly like working correctly.
            if note and note.lower() not in {"nan", "none", "false", "0", "no", "-"}:
                reasons.append(f"needs_confirmation={note[:120]!r}")
        if reasons:
            bad.append((idx, "; ".join(reasons)))
    return bad


def _standard_exposure(df, exposure_cols):
    """The cohort's modal exposure triple, as a tuple. None if it cannot be found."""
    triples = []
    for _, row in df.iterrows():
        values = tuple(row[c] for c in exposure_cols)
        if any(_is_missing(v) for v in values):
            continue
        triples.append(tuple(round(float(v), 4) for v in values))
    if not triples:
        return None
    return statistics.mode(triples)


def _exposure_is_standard(df, exposure_cols, standard):
    out = []
    for _, row in df.iterrows():
        values = [row[c] for c in exposure_cols]
        if any(_is_missing(v) for v in values):
            out.append(None)
            continue
        out.append(all(
            math.isclose(float(v), s, rel_tol=EXPOSURE_RTOL)
            for v, s in zip(values, standard)))
    return out


def _has_duplicate_physical_sections(df, tube_col):
    """True if any (tube, physical section) appears more than once in *df*.

    Uses the PHYSICAL section identity where the manifest supplies it, falling back
    to the scanner's label. The distinction matters: tube 51's rescan labels its two
    sections 01/02 while they are physically the original slide's 03/04, so the
    scanner label alone both misses real duplicates and invents false ones.
    """
    label_col = _column(df, "physical_section_label") or _column(df, "section_label")
    if tube_col is None or label_col is None:
        return False
    return bool(df.duplicated(subset=[tube_col, label_col]).any())


def split_manifest(df, codes: Mapping, *, extra_blinded_columns: Sequence[str] = (),
                   collapse_duplicate_scans: bool = True):
    """Split one manifest into (private provenance, blinded analysis).

    The private frame keeps everything and gains a ``code`` column. The blinded frame
    keeps only the columns named in :data:`BLINDED_COLUMNS` (plus
    ``extra_blinded_columns``), never a column whose name matches
    :data:`NEVER_BLINDED`, and gains ``exposure_is_standard`` in place of the raw
    exposure times. Everything else is dropped and listed in
    ``blinded.attrs["blinding_split_report"]["dropped_columns"]``.

    Two things are done that are easy to forget and both leak on their own:

    * **Row order is destroyed.** The blinded frame is sorted by ``(code,
      section_label)`` and re-indexed. A manifest left in file-iteration order is in
      tube order, which is group order, which means the row numbers alone would
      rebuild the design.
    * **Pixel size is rounded** to :data:`PIXEL_SIZE_DECIMALS`. At full precision it
      varies in the fifth decimal per series and is therefore a fingerprint that joins
      straight back to ``work/meta/series_metadata.csv``, which carries ``tube_id``.

    Args:
        df: The analysis manifest, one row per section (or per section and region).
            Must carry ``tube_id``.
        codes: ``{tube_id: code}`` from :func:`generate_codes`, covering every tube in
            ``df``.
        extra_blinded_columns: Additional column names to let through, for a caller
            who has thought about it. Names matching :data:`NEVER_BLINDED` are still
            refused.
        collapse_duplicate_scans: Keep only the rows the manifest marks
            ``scan_is_preferred``. **On by default**, because a blinded analysis
            manifest containing two rows for one physical section is not a usable
            artefact: anything addressing a section by ``(code, section)`` gets two
            answers, and in this cohort the two rows can carry opposite conditions.
            It is also a tell -- only two animals were re-acquired, so the presence
            of a duplicate identifies them. The repeat scans remain in the PRIVATE
            manifest, which is where the ADR-0012 repeatability check reads them.
            Pass ``False`` (``./ihc blind --keep-all-scans``) only if you have a
            specific reason and have thought about both consequences.

    Returns:
        ``(private_df, blinded_df)``.

    Raises:
        UnresolvedConditionError: any row's condition is unresolved, or any row
            carries a non-empty ``needs_confirmation``. Tube 37 is the live case.
        KeyError: ``df`` has no ``tube_id`` column, or a tube in ``df`` has no code.
        ValueError: ``extra_blinded_columns`` names a forbidden column.
    """
    pd = _pandas()
    tube_col = _column(df, "tube_id")
    if tube_col is None:
        raise KeyError(
            "manifest has no 'tube_id' column, so it cannot be coded. Columns: "
            f"{list(df.columns)}")

    collapsed_note = None
    if collapse_duplicate_scans:
        preferred_col = _column(df, "scan_is_preferred", "is_preferred_scan", "preferred")
        if preferred_col is None:
            # No such column. That is only a problem if there is actually something to
            # collapse -- a manifest with one scan per section (or any synthetic frame)
            # needs no rule, so silently do nothing. Raise only when duplicates exist
            # and there is no basis for choosing between them, because picking one
            # arbitrarily is exactly the kind of quiet choice that produces a wrong
            # number nobody can trace.
            if _has_duplicate_physical_sections(df, tube_col):
                raise ValueError(
                    "collapse_duplicate_scans=True and this manifest has more than one "
                    "row for the same physical section, but no 'scan_is_preferred' "
                    "column to say which scan to keep. That choice belongs to the "
                    "manifest, not to this function. Rebuild the manifest, or pass "
                    "collapse_duplicate_scans=False and resolve it downstream.")
            collapse_duplicate_scans = False
    # The collapse applies to the BLINDED frame only. The private provenance manifest
    # is the archival record of what was imaged -- it must keep every scan of every
    # section, or the repeat scans that ADR-0012 uses for the imaging-repeatability
    # check vanish from the one file the custodian is meant to freeze and back up.
    # Earlier this filtered both frames together and silently dropped the original
    # scans of tubes 51 and 60 from the provenance record.
    blinded_keep = None
    if collapse_duplicate_scans:
        blinded_keep = df[preferred_col].map(lambda v: str(v).strip().lower() in
                                             {"true", "1", "yes", "y"})
        collapsed_note = (
            f"collapse_duplicate_scans: the blinded manifest keeps "
            f"{int(blinded_keep.sum())} of {len(df)} row(s) marked {preferred_col}. "
            "The rest are repeat scans of the same physical section; they remain in "
            "the private provenance manifest and in manifest.csv.")

    unresolved = _unresolved_rows(df)
    if unresolved:
        detail = []
        for idx, reason in unresolved[:10]:
            tube = df.at[idx, tube_col]
            detail.append(f"    row {idx} (tube {tube}): {reason}")
        more = "" if len(unresolved) <= 10 else f"\n    ... and {len(unresolved) - 10} more"
        raise UnresolvedConditionError(
            f"{len(unresolved)} manifest row(s) have an unresolved condition:\n"
            + "\n".join(detail) + more
            + "\n  Blinding these would hide the contradiction under a code and the "
              "numbers would still look plausible. Exclude them from the analysis "
              "manifest and report them; resolve at the bench, not from the pixels "
              "(ADR-0003).")

    missing = sorted({int(t) for t in df[tube_col] if not _is_missing(t)}
                     - {int(k) for k in codes})
    if missing:
        raise KeyError(f"no code for tube(s) {missing}; regenerate codes for the full cohort")

    code_by_tube = {int(k): v for k, v in codes.items()}

    # ---- private: everything, plus the code -------------------------------------
    private = df.copy()
    private.insert(0, "code", [code_by_tube[int(t)] for t in private[tube_col]])

    # ---- blinded: an allow-list, and everything else is reported and dropped -----
    allow = list(BLINDED_COLUMNS)
    for extra in extra_blinded_columns:
        if _forbidden_name(extra):
            raise ValueError(
                f"{extra!r} matches NEVER_BLINDED ({_forbidden_name(extra)!r}) and "
                "cannot be added to the blinded manifest")
        if extra not in allow:
            allow.append(extra)

    blinded = pd.DataFrame(index=df.index)
    blinded["code"] = private["code"].values

    kept, dropped, notes = ["code"], [], ([collapsed_note] if collapsed_note else [])
    for column in df.columns:
        if column == tube_col:
            dropped.append(column)
            continue
        forbidden = _forbidden_name(column)
        if forbidden:
            dropped.append(column)
            continue
        if column in allow:
            values = df[column]
            if str(column).lower() == "pixel_size_um":
                values = values.astype(float).round(PIXEL_SIZE_DECIMALS)
                notes.append(f"pixel_size_um rounded to {PIXEL_SIZE_DECIMALS} decimals "
                             "(full precision fingerprints the series)")
            blinded[column] = values.values
            kept.append(column)
        else:
            dropped.append(column)

    # ---- exposure -> a single boolean -------------------------------------------
    exposure_cols = [c for c in df.columns
                     if re.search(r"exposure", str(c), re.IGNORECASE)
                     and re.search(r"dapi|fitc|cy3|_ms$|_us$", str(c), re.IGNORECASE)]
    if "exposure_is_standard" not in blinded.columns and exposure_cols:
        standard = _standard_exposure(df, exposure_cols)
        if standard is not None:
            blinded["exposure_is_standard"] = _exposure_is_standard(
                df, exposure_cols, standard)
            kept.append("exposure_is_standard")
            notes.append(
                f"exposure_is_standard derived from the cohort modal triple "
                f"{dict(zip([str(c) for c in exposure_cols], standard))}; raw exposure "
                "columns dropped. NOTE this boolean is False for exactly the two "
                "animals the spec names in writing -- it coarsens the leak, it does "
                "not remove it.")

    # ---- make every row addressable, without saying which scan it is -------------
    #
    # Two scans of one physical slide give two rows with the same
    # (code, section_label) and -- because the rescan covers the positive box only --
    # OPPOSITE conditions. Left alone that is a correctness hazard of exactly the kind
    # spec §2 is about: a reader cannot tell which row is the negative control.
    #
    # The obvious fix, a `scan`/`is_rescan` column, is not available. Only two animals
    # were re-acquired; ADR-0006 and ADR-0012 name them in writing; so a column marking
    # two of thirty animals as rescans resolves to those tubes for any reader who has
    # the project documents -- which is everyone, including the person who must stay
    # blinded. "Rescan status identifies no animal on its own" is true only for a reader
    # with no documents, and that is not the threat model here.
    #
    # So: a neutral within-section replicate index, ordered by row CONTENT rather than
    # by scan order, so it says "there are two of these" without saying which came from
    # which session. The animals with duplicates are still identifiable -- duplication
    # is itself the tell -- which is why the audit raises it and why
    # collapse_duplicate_scans exists.
    label_col = "section_label" if "section_label" in blinded.columns else None
    if (label_col is not None and "replicate" not in blinded.columns
            and blinded.duplicated(["code", label_col]).any()):
        rest = blinded.drop(columns=["code", label_col])
        if rest.columns.empty:
            # Nothing to order by. Fall back to position: with no other column the two
            # rows are identical anyway, so the index carries no information either way.
            order = blinded.groupby(["code", label_col]).cumcount() + 1
        else:
            content = rest.astype(str).agg("|".join, axis=1)
            order = content.groupby([blinded["code"], blinded[label_col]]).rank(method="first")
        blinded["replicate"] = order.astype(int)
        kept.append("replicate")
        notes.append(
            "replicate added: some sections appear more than once (repeat scans of the "
            "same physical slide). The index is ordered by row content, not by scan "
            "order, so it does not say which acquisition is which -- but the presence "
            "of duplicates identifies those animals by itself. Prefer "
            "collapse_duplicate_scans=True (./ihc blind --preferred-scan-only).")

    # ---- collapse repeat scans, in the BLINDED frame only ------------------------
    if blinded_keep is not None:
        blinded = blinded[blinded_keep.values]

    # ---- destroy row order ------------------------------------------------------
    sort_keys = [c for c in ("code", "section_label", "replicate", "region",
                             "series_index") if c in blinded.columns]
    blinded = blinded.sort_values(sort_keys, kind="mergesort").reset_index(drop=True)
    private = private.reset_index(drop=True)

    blinded.attrs["blinding_split_report"] = {
        "kept_columns": kept,
        "dropped_columns": dropped,
        "notes": notes,
        "sorted_by": sort_keys,
    }
    return private, blinded


def _forbidden_name(column) -> str | None:
    """The NEVER_BLINDED token that `column` matches, or None."""
    name = str(column).lower()
    for token in NEVER_BLINDED:
        if token in name:
            return token
    return None


# --------------------------------------------------------------------------- #
# Small statistics, stdlib only
# --------------------------------------------------------------------------- #


def _ranks(values: Sequence[float]) -> list[float]:
    """Average ranks, 1-based."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = mean_rank
        i = j + 1
    return out


def _pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return 0.0
    return sxy / math.sqrt(sxx * syy)


def _spearman(xs: Sequence[float], ys: Sequence[float]) -> float:
    return _pearson(_ranks(list(xs)), _ranks(list(ys)))


def _adjacent_same_group(groups: Sequence[str]) -> int:
    return sum(1 for a, b in zip(groups, groups[1:]) if a == b)


def _mc_pvalues(tubes: Sequence[int], code_ranks: Sequence[float],
                groups: Sequence[str], n_perm: int = AUDIT_MC_PERMUTATIONS):
    """Monte Carlo nulls for both code-order tests, from one shared permutation stream.

    Returns ``(rho, p_rho_two_sided, adjacency, p_adj_upper, p_adj_lower, null_mean)``.

    The null is "the code order is a uniformly random permutation of the animals",
    which is precisely what :func:`generate_codes` claims to produce. The adjacency
    statistic -- how many neighbouring animals in code order share a treatment group --
    is the test with real content: a coding scheme that keeps the treatment blocks
    contiguous scores near ``n - 1`` where chance scores about 7.7 for this cohort,
    and it fails here even in the cases where the rank correlation happens to look
    innocent (a block-reversed or block-rotated order, for instance).
    """
    rho = _spearman(tubes, code_ranks)
    obs_adj = _adjacent_same_group([g for _, g in
                                    sorted(zip(code_ranks, groups), key=lambda t: t[0])])

    rng = random.Random(AUDIT_MC_SEED)
    order = list(range(len(tubes)))
    ge_rho = ge_adj = le_adj = 0
    null_adj_total = 0
    tube_ranks = _ranks(list(tubes))
    for _ in range(n_perm):
        rng.shuffle(order)
        permuted_ranks = [code_ranks[i] for i in order]
        if abs(_pearson(tube_ranks, _ranks(permuted_ranks))) >= abs(rho) - 1e-12:
            ge_rho += 1
        adj = _adjacent_same_group([groups[i] for i in
                                    sorted(range(len(order)), key=lambda k: permuted_ranks[k])])
        null_adj_total += adj
        ge_adj += adj >= obs_adj
        le_adj += adj <= obs_adj
    return (rho,
            (1 + ge_rho) / (n_perm + 1),
            obs_adj,
            (1 + ge_adj) / (n_perm + 1),
            (1 + le_adj) / (n_perm + 1),
            null_adj_total / n_perm)


def _partition(labels_by_animal: Mapping) -> frozenset:
    cells = {}
    for animal, label in labels_by_animal.items():
        cells.setdefault(label, set()).add(animal)
    return frozenset(frozenset(c) for c in cells.values())


# --------------------------------------------------------------------------- #
# The audit
# --------------------------------------------------------------------------- #


def audit_blinded(blinded_df, codes: Mapping, private_df) -> list:
    """Try to break the blinding. Return every leak found, worst first.

    This is the function that matters. Everything else in the module is bookkeeping;
    this is the part that decides whether the bookkeeping worked. It has access to the
    truth (``private_df``) precisely so that it can attack the blinded artefact the way
    someone holding the bench record would.

    What it runs:

    1. **Tube IDs.** Every column and every column *name*, searched for each tube ID as
       a value or as a digit-bounded substring. Digit-bounded so that ``15360`` does
       not read as tube 36.
    2. **Group and arm labels.** Full labels and their distinctive tokens, plus a
       generic treatment vocabulary, case-insensitively. The phrases "negative
       control" and "positive control" are masked out first, so ``has_negative_control``
       is not mistaken for the group ``Control IP (vehicle)``.
    3. **Paths and file names.** Absolute and relative paths, UNC and drive-letter
       paths, known image and project extensions, and this cohort's ``Image_NN`` /
       ``_Image_NN_`` stems.
    4. **Is the code order order-preserving?** Spearman between tube ID and code rank,
       with a Monte Carlo p-value; **and** a permutation test on the number of adjacent
       same-group pairs in code order, which is the test with real content. Also the
       row order of the blinded frame itself, which is an ordering nobody thinks about
       and which is in tube order unless it was deliberately destroyed.
    5. **Does any column recover the design?** For every column, the partition of
       animals induced by its values is compared with the true group partition:

       * exact recovery, or a coarsening in which every cell is group-pure -> ``LEAK``
       * monotone in tube ID (|rho| >= 0.9) -> ``LEAK`` -- this is the acquisition
         timestamp case, and note that ISO timestamps sort correctly as strings
       * one-to-one with animal -> ``LEAK``: it need not encode group itself, because
         it joins straight back to any un-blinded table that carries ``tube_id``
       * a value held by only one or two animals -> ``RISK``: it does not name the
         animal but it singles it out, which is enough for anyone with the bench record.
         This is the raw-exposure case, and also the case for the ``exposure_is_standard``
         boolean that replaces it.

    Note what test 5's *first* clause does **not** catch on its own, since the spec
    expected it to: raw exposure splits the cohort 29/1/1, which is not the group
    partition, and acquisition timestamps are all distinct, which is not the group
    partition either. Exact partition equality alone would pass both. They are caught
    by the monotonicity, uniqueness and rare-value clauses instead -- which is why the
    test is a family rather than a single comparison.

    Args:
        blinded_df: the frame that is about to be handed out.
        codes: ``{tube_id: code}``.
        private_df: the custodian's frame, carrying ``tube_id``, ``group`` and ``code``.

    Returns:
        A list of strings, each starting ``LEAK:`` or ``RISK:``, ordered leaks first.
        Empty means nothing was found -- which is not the same as "there is nothing".

    Side effect:
        The measured statistics (rho, its p-value, the adjacency count and its null)
        are written to ``blinded_df.attrs["blinding_audit_stats"]``, because a passing
        test still has a number worth recording. They are kept out of the return value
        so that "empty list" keeps meaning "no finding".
    """
    leaks, risks = [], []
    stats = {}
    code_col = _column(blinded_df, "code", "coded_id", "code_id")
    if code_col is None:
        raise KeyError("blinded manifest has no 'code' column")

    tube_by_code = {v: int(k) for k, v in codes.items()}
    tubes_all = sorted(tube_by_code.values())

    p_tube = _column(private_df, "tube_id")
    p_group = _column(private_df, "group")
    p_code = _column(private_df, "code")
    if p_tube is None or p_code is None:
        raise KeyError("private manifest needs 'tube_id' and 'code' columns")

    group_by_tube = {}
    if p_group is not None:
        for _, row in private_df.iterrows():
            group_by_tube[int(row[p_tube])] = _text(row[p_group]).strip()

    # ---- 0. structural ----------------------------------------------------------
    blinded_codes = {c for c in blinded_df[code_col] if not _is_missing(c)}
    unknown = sorted(blinded_codes - set(tube_by_code))
    if unknown:
        leaks.append(f"LEAK: blinded manifest carries code(s) with no key entry: {unknown}")
    if len(set(codes.values())) != len(codes):
        leaks.append("LEAK: the code map is not one-to-one -- two animals share a code")
    for column in blinded_df.columns:
        forbidden = _forbidden_name(column)
        if forbidden:
            leaks.append(
                f"LEAK: column {column!r} matches the never-blinded token "
                f"{forbidden!r}. Every entry in NEVER_BLINDED is a leak vector verified "
                "present in this dataset, so the name alone is disqualifying -- the "
                "values do not get a hearing.")

    # ---- 1. tube IDs ------------------------------------------------------------
    for column in blinded_df.columns:
        hits = _tube_hits_in_name(str(column), tubes_all)
        if hits:
            leaks.append(f"LEAK: column name {column!r} contains tube ID(s) {hits}")
        hits, severity = _tube_hits_in_values(blinded_df[column], tubes_all)
        if hits and severity == "LEAK":
            leaks.append(
                f"LEAK: column {column!r} contains tube ID(s) {hits} as values")
        elif hits:
            risks.append(
                f"RISK: column {column!r} contains value(s) {hits} that coincide with "
                "tube ID(s); probably a coincidence, check it is")

    # ---- 2. group and arm labels ------------------------------------------------
    vocabulary = _group_vocabulary(private_df, p_group)
    for column in blinded_df.columns:
        found = _vocabulary_hits(str(column), vocabulary)
        if found:
            leaks.append(f"LEAK: column name {column!r} contains group/arm term(s) {found}")
        found = set()
        for value in blinded_df[column]:
            if isinstance(value, str):
                found |= set(_vocabulary_hits(value, vocabulary))
        if found:
            leaks.append(
                f"LEAK: column {column!r} contains group/arm term(s) {sorted(found)} in its values")

    # ---- 3. paths and file names ------------------------------------------------
    for column in blinded_df.columns:
        if _PATHISH.search(str(column)) or _forbidden_name(column) in ("path", "file", "dir", "folder"):
            leaks.append(f"LEAK: column name {column!r} looks like a path or file reference")
        examples = [v for v in blinded_df[column]
                    if isinstance(v, str) and _PATHISH.search(v)]
        if examples:
            leaks.append(
                f"LEAK: column {column!r} contains path/filename-like value(s), "
                f"e.g. {examples[0][:80]!r} ({len(examples)} row(s))")

    # ---- 4. is the code order order-preserving? ---------------------------------
    coded_tubes = sorted(tube_by_code.values())
    ranks_by_code = {code: rank for rank, code in enumerate(sorted(codes.values()))}
    code_ranks = [ranks_by_code[codes[t]] for t in coded_tubes]
    groups_seq = [group_by_tube.get(t, "?") for t in coded_tubes]

    if len(coded_tubes) >= 5:
        rho, p_rho, adj, p_adj_hi, p_adj_lo, null_mean = _mc_pvalues(
            coded_tubes, code_ranks, groups_seq)
        stats.update({
            "n_animals": len(coded_tubes),
            "spearman_tube_vs_code_rank": round(rho, 4),
            "spearman_mc_p_two_sided": round(p_rho, 5),
            "adjacent_same_group_pairs": adj,
            "adjacent_same_group_null_mean": round(null_mean, 3),
            "adjacent_mc_p_upper": round(p_adj_hi, 5),
            "adjacent_mc_p_lower": round(p_adj_lo, 5),
            "mc_permutations": AUDIT_MC_PERMUTATIONS,
            "mc_seed": AUDIT_MC_SEED,
        })
        if abs(rho) >= MAX_ABS_SPEARMAN or p_rho <= MC_ALPHA:
            leaks.append(
                f"LEAK: code order is order-preserving in tube ID: Spearman rho="
                f"{rho:+.3f} (|rho| limit {MAX_ABS_SPEARMAN}, Monte Carlo two-sided "
                f"p={p_rho:.4f}). Tube IDs run in contiguous treatment blocks, so this "
                "reproduces the group structure.")
        elif p_rho <= MC_ADVISORY_ALPHA:
            risks.append(
                f"RISK: the draw is more monotone in tube ID than most: Spearman "
                f"rho={rho:+.3f}, Monte Carlo two-sided p={p_rho:.3f}. Below the "
                f"failure threshold, but a seed is free -- re-drawing until the order "
                "test is unremarkable is restricted randomisation, not cheating. "
                "Re-draw BEFORE anything is coded, and record which seed was used.")
        if group_by_tube:
            if p_adj_hi <= MC_ALPHA:
                leaks.append(
                    f"LEAK: treatment blocks stay contiguous in code order: "
                    f"{adj} adjacent same-group pairs against a null mean of "
                    f"{null_mean:.2f} (Monte Carlo upper-tail p={p_adj_hi:.4f}). "
                    "A code sequence that keeps the blocks together has failed.")
            elif p_adj_lo <= MC_ALPHA:
                risks.append(
                    f"RISK: code order is *anti*-blocked -- {adj} adjacent same-group "
                    f"pairs against a null mean of {null_mean:.2f} (lower-tail "
                    f"p={p_adj_lo:.4f}). Structured rather than random; check the seed "
                    "was not chosen to produce it.")

    # row order of the blinded file itself
    row_tubes = [tube_by_code.get(c) for c in blinded_df[code_col]]
    row_tubes = [t for t in row_tubes if t is not None]
    if len(set(row_tubes)) >= 5:
        rho_rows = _spearman(list(range(len(row_tubes))), row_tubes)
        stats["spearman_row_order_vs_tube"] = round(rho_rows, 4)
        if abs(rho_rows) >= MAX_ABS_SPEARMAN:
            leaks.append(
                f"LEAK: the ROW ORDER of the blinded manifest is monotone in tube ID "
                f"(Spearman rho={rho_rows:+.3f}). Row number alone rebuilds the design; "
                "sort by code before writing.")

    # ---- 5. does any column recover the design? ---------------------------------
    for column in blinded_df.columns:
        if column == code_col:
            continue
        per_animal = {}
        for code, value in zip(blinded_df[code_col], blinded_df[column]):
            tube = tube_by_code.get(code)
            if tube is None:
                continue
            per_animal.setdefault(tube, set()).add(_text(value).strip())
        signature = {t: "|".join(sorted(v)) for t, v in per_animal.items()}
        distinct = set(signature.values())
        if len(distinct) <= 1:
            continue                                    # constant: leaks nothing

        n_animals = len(signature)
        # 5e -- one-to-one with animal
        if len(distinct) == n_animals and n_animals >= 3:
            leaks.append(
                f"LEAK: column {column!r} takes a distinct value for every animal. It "
                "need not encode group: it is a join key back to any un-blinded table "
                "that carries tube_id (work/meta/series_metadata.csv does).")
            continue

        if group_by_tube:
            group_partition = _partition({t: group_by_tube.get(t, "?") for t in signature})
            column_partition = _partition(signature)
            if column_partition == group_partition:
                leaks.append(
                    f"LEAK: column {column!r} exactly recovers the treatment-group "
                    "partition -- grouping animals by its value reproduces the design.")
                continue
            pure = all(len({group_by_tube.get(t, "?") for t in cell}) == 1
                       for cell in column_partition)
            if pure and len(column_partition) < n_animals:
                leaks.append(
                    f"LEAK: column {column!r} never mixes treatment groups "
                    f"({len(column_partition)} value-cells over {n_animals} animals, "
                    "every cell group-pure). It is a refinement of the design, which "
                    "is as good as the design.")
                continue

            # 5c -- monotone in tube ID
            rho_col = _column_rho(signature)
            if rho_col is not None and abs(rho_col) >= COLUMN_MONOTONE_RHO:
                leaks.append(
                    f"LEAK: column {column!r} is monotone in tube ID (Spearman "
                    f"rho={rho_col:+.3f}). Scan order is tube order is group order, so "
                    "an ordering is a group label. This is the acquisition-timestamp case.")
                continue

        # 5d -- rare values single out individual animals
        counts = {}
        for tube, value in signature.items():
            counts.setdefault(value, []).append(tube)
        rare = {v: len(ts) for v, ts in counts.items()
                if len(ts) <= RARE_VALUE_MAX_ANIMALS}
        if rare and n_animals >= 8:
            risks.append(
                f"RISK: column {column!r} has value(s) {sorted(rare)} carried by only "
                f"{sorted(rare.values())} of the {n_animals} animals. That does not "
                "name an animal, but whoever holds the bench record can: a value almost "
                "nobody shares is a pointer, and the exceptional slides in this cohort "
                "are named in writing in the spec and the ADRs. Either drop the column "
                "or accept it and record why.")

    # ---- 5e2. duplicate coded rows: the rescan tell ------------------------------
    label_col = _column(blinded_df, "section_label")
    if label_col is not None:
        seen, duplicated = {}, set()
        for code, label in zip(blinded_df[code_col], blinded_df[label_col]):
            key = (code, _text(label))
            if key in seen:
                duplicated.add(code)
            seen[key] = True
        has_replicate = _column(blinded_df, "replicate") is not None
        if duplicated and len(duplicated) <= RARE_VALUE_MAX_ANIMALS:
            risks.append(
                f"RISK: {len(duplicated)} animal(s) appear more than once for the same "
                "section_label. In this cohort that means a repeat scan; the "
                "re-acquired slides are named in writing in ADR-0006 and ADR-0012; so "
                "the duplication itself points at those animals no matter what the "
                "extra rows are called. "
                + ("A neutral 'replicate' index makes them addressable, but it does not "
                   "make them anonymous. " if has_replicate else
                   "There is also no discriminator, so a section cannot be addressed by "
                   "(code, section_label) at all -- and the two rows can carry OPPOSITE "
                   "conditions, which is the failure mode spec §2 exists to prevent. ")
                + "Blind the preferred scan only (./ihc blind --preferred-scan-only) and "
                  "keep the repeat-scan repeatability check on the custodian side.")
        elif duplicated:
            risks.append(
                f"RISK: {len(duplicated)} animals have more than one row per "
                "section_label"
                + ("." if has_replicate else
                   ", so a section cannot be addressed unambiguously by "
                   "(code, section_label). Downstream needs a rule."))

    # ---- 5f. the JOINT signature, which can identify where no single column does --
    joint = {}
    payload_cols = [c for c in blinded_df.columns if c != code_col]
    if payload_cols:
        for code, *values in zip(blinded_df[code_col],
                                 *(blinded_df[c] for c in payload_cols)):
            tube = tube_by_code.get(code)
            if tube is None:
                continue
            joint.setdefault(tube, set()).add("|".join(_text(v).strip() for v in values))
        signature = {t: "||".join(sorted(v)) for t, v in joint.items()}
        counts = {}
        for tube, value in signature.items():
            counts.setdefault(value, []).append(tube)
        singled = sum(len(ts) for ts in counts.values()
                      if len(ts) <= RARE_VALUE_MAX_ANIMALS)
        if len(counts) == len(signature) and len(signature) >= 3:
            leaks.append(
                "LEAK: the whole row-set is a fingerprint -- every animal has a unique "
                "combination of blinded values, so the manifest joins one-to-one against "
                "any un-blinded table even though no single column does.")
        elif singled:
            risks.append(
                f"RISK: {singled} of {len(signature)} animals have a combination of "
                "blinded values shared with no more than "
                f"{RARE_VALUE_MAX_ANIMALS - 1} other animal(s). Individually every "
                "column is safe; jointly they narrow. Some of this is unavoidable -- "
                "three-section slides and 4+0 slides really are different from the rest "
                "-- so it is recorded rather than fixed. Check it is not more than the "
                "cohort structure explains.")
        stats["joint_signature_cells"] = len(counts)
        stats["joint_signature_singled_animals"] = singled

    try:
        blinded_df.attrs["blinding_audit_stats"] = stats
    except Exception:                                          # noqa: BLE001
        pass                                   # a frame that will not take attrs is fine
    return leaks + risks


def _column_rho(signature: Mapping) -> float | None:
    """Spearman between tube ID and a per-animal column value, or None if unrankable.

    Numeric values are used directly; anything else is ranked lexicographically, which
    is what makes ISO-8601 timestamps -- the most dangerous case -- fall out correctly.
    """
    tubes = sorted(signature)
    values = [signature[t] for t in tubes]
    try:
        numeric = [float(v) for v in values]
    except (TypeError, ValueError):
        numeric = None
    if numeric is None:
        order = {v: i for i, v in enumerate(sorted(set(values)))}
        numeric = [order[v] for v in values]
    if len(set(numeric)) < 3:
        return None
    return _spearman(tubes, numeric)


def _tube_hits_in_name(name: str, tubes: Sequence[int]) -> list:
    masked = _NUMERIC_NOISE.sub(" ", str(name))
    return [t for t in tubes if re.search(rf"(?<!\d){t}(?!\d)", masked)]


def _tube_hits_in_values(series, tubes: Sequence[int]):
    """(hits, severity) for one column. Strings are searched, numbers compared."""
    tube_set = set(tubes)
    hits = set()
    is_float = False
    for value in series:
        if _is_missing(value):
            continue
        if isinstance(value, str):
            masked = _NUMERIC_NOISE.sub(" ", value)
            hits |= {t for t in tubes if re.search(rf"(?<!\d){t}(?!\d)", masked)}
        elif isinstance(value, bool):
            continue
        elif isinstance(value, (int,)) or (hasattr(value, "is_integer") and float(value).is_integer()):
            is_float = is_float or isinstance(value, float)
            if int(value) in tube_set:
                hits.add(int(value))
    if not hits:
        return [], None
    # A float column that happens to land on one tube number is probably arithmetic;
    # two or more, or any string hit, is not.
    severity = "RISK" if (is_float and len(hits) == 1) else "LEAK"
    return sorted(hits), severity


def _group_vocabulary(private_df, group_col) -> list:
    """Full group/arm labels plus their distinctive tokens, lower-cased."""
    vocabulary = set(_TREATMENT_WORDS)
    columns = [c for c in (group_col, _column(private_df, "arm")) if c is not None]
    for column in columns:
        for value in private_df[column]:
            text = _text(value).strip().lower()
            if not text:
                continue
            vocabulary.add(text)
            for token in re.findall(r"[a-z]{3,}", text):
                # "control" and "extra" are the vocabulary of negative controls and of
                # ordinary English; they are handled by masking benign phrases instead
                # of by dropping the token, so a bare "Control" value is still caught.
                vocabulary.add(token)
    return sorted(vocabulary)


def _vocabulary_hits(text: str, vocabulary: Sequence[str]) -> list:
    masked = _BENIGN_CONTROL.sub(" ", str(text)).lower()
    masked = masked.replace("_", " ")
    return sorted({term for term in vocabulary
                   if re.search(rf"(?<![a-z]){re.escape(term)}(?![a-z])", masked)})


# --------------------------------------------------------------------------- #
# Writing
# --------------------------------------------------------------------------- #


def _check_custodian_dir(path: Path) -> None:
    """Refuse a custodian directory that git could reach.

    Cloud sync is deliberately allowed. An earlier version refused any path under
    Dropbox on the grounds that it might be shared; the PI's decision is that the
    project folder is not shared with the delineator, and that syncing the key is
    how it gets backed up. A key that exists in exactly one place on one laptop is
    the bigger risk: lose it and every measurement made under it is orphaned.

    Git is still refused, because a pushed key is public and cannot be recalled.
    """
    resolved = path.expanduser().resolve()
    for ancestor in [resolved] + list(resolved.parents):
        if (ancestor / ".git").exists():
            raise CustodianPathError(
                f"custodian directory is inside a git repository ({ancestor}): "
                f"{resolved}\n  The key would be committable and pushable. "
                "ADR-0001: the custodian tree lives outside git.")


def _atomic_write(path: Path, text: str, mode: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            fh.write(text)
        os.chmod(tmp, mode)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def write_blinded(private_df, blinded_df, codes: Mapping, seed, *,
                  custodian_dir, work_dir, force: bool = False,
                  audit: Sequence[str] | None = None) -> dict:
    """Write the key and both manifests, and refuse if the audit finds a hard leak.

    Three artefacts go to ``custodian_dir`` (created mode 0700, files 0600): the key
    CSV, a key JSON carrying **the seed** and enough provenance to regenerate the
    mapping, and the private provenance manifest with the seed and a UTC timestamp
    stamped onto every row. One artefact plus a report goes to ``work_dir``: the
    blinded analysis manifest, and a JSON report that lists the audit findings and the
    dropped columns but contains **no seed and no code mapping**.

    The audit runs here whether or not the caller ran it. Any ``LEAK:`` finding stops
    the blinded manifest being written at all, and ``force`` does not override that --
    ``force`` overrides only the refusal to overwrite an existing key.

    Args:
        private_df: from :func:`split_manifest`.
        blinded_df: from :func:`split_manifest`.
        codes: ``{tube_id: code}``.
        seed: the seed that produced ``codes``. Recorded with the key, never in
            ``work_dir``.
        custodian_dir: mode-0700 directory outside git and outside Dropbox.
        work_dir: where the blinded manifest goes; may be shared storage.
        force: overwrite an existing key. Without it, an existing key is an error --
            silently rewriting a key orphans every coded artefact already produced
            under the old one.
        audit: a previously computed audit, to avoid running it twice.

    Returns:
        A dict of paths written, counts, dropped columns, and the audit findings.

    Raises:
        FileExistsError: a key exists and ``force`` is false.
        BlindingLeakError: the audit found a ``LEAK:``.
        CustodianPathError: ``custodian_dir`` is inside git or inside Dropbox.
    """
    custodian = Path(custodian_dir).expanduser()
    work = Path(work_dir).expanduser()
    _check_custodian_dir(custodian)
    if custodian.resolve() == work.expanduser().resolve():
        raise CustodianPathError(
            "custodian_dir and work_dir are the same directory; the key would sit "
            "beside the blinded manifest")

    # The seed about to be recorded must actually regenerate the mapping about to be
    # recorded. If it does not, the key is a note of a number that opens nothing, and
    # nobody would find out until the day the cohort had to be un-blinded.
    seed = normalise_seed(seed)
    regenerated = generate_codes(sorted(int(t) for t in codes), seed=seed)
    if {int(k): v for k, v in codes.items()} != regenerated:
        raise ValueError(
            "the seed does not regenerate the code mapping. Either the codes were not "
            "produced by generate_codes(tube_ids, seed=seed) over exactly this set of "
            "animals, or the seed passed here is not the one that produced them. "
            "Recording it would leave a key that opens nothing.")

    findings = list(audit) if audit is not None else audit_blinded(blinded_df, codes, private_df)
    hard = [f for f in findings if f.startswith("LEAK:")]
    soft = [f for f in findings if not f.startswith("LEAK:")]

    key_csv = custodian / _KEY_CSV_NAME
    key_json = custodian / _KEY_JSON_NAME
    provenance_csv = custodian / _PROVENANCE_CSV_NAME
    blinded_csv = work / _BLINDED_CSV_NAME
    report_json = work / _REPORT_JSON_NAME

    existing = [p for p in (key_csv, key_json) if p.exists()]
    if existing and not force:
        raise FileExistsError(
            f"a blinding key already exists: {', '.join(str(p) for p in existing)}\n"
            "  Overwriting it orphans every coded artefact produced under the old key "
            "-- QuPath projects, annotations, measurements. Pass force=True (./ihc "
            "blind --force) only if you are certain nothing downstream has been made yet.")

    if hard:
        raise BlindingLeakError(
            f"the audit found {len(hard)} hard leak(s); nothing was written.\n  "
            + "\n  ".join(hard)
            + "\n  Fix the manifest columns and re-run. force does not override this.")

    custodian.mkdir(parents=True, exist_ok=True)
    os.chmod(custodian, 0o700)
    work.mkdir(parents=True, exist_ok=True)

    created = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # ---- private provenance -----------------------------------------------------
    seed = normalise_seed(seed)
    provenance = private_df.copy()
    provenance["blinding_seed"] = str(seed)
    provenance["blinded_at_utc"] = created
    _atomic_write(provenance_csv, provenance.to_csv(index=False), 0o600)

    # ---- the key ----------------------------------------------------------------
    pd = _pandas()
    key_columns = {"tube_id": [int(t) for t in sorted(codes)],
                   "code": [codes[t] for t in sorted(codes)]}
    for name in ("group", "arm"):
        column = _column(private_df, name)
        if column is None:
            continue
        by_tube = {}
        tube_col = _column(private_df, "tube_id")
        for _, row in private_df.iterrows():
            by_tube[int(row[tube_col])] = row[column]
        key_columns[name] = [by_tube.get(t) for t in sorted(codes)]
    key_frame = pd.DataFrame(key_columns)
    _atomic_write(key_csv, key_frame.to_csv(index=False), 0o600)

    split_report = dict(blinded_df.attrs.get("blinding_split_report", {}))
    audit_stats = dict(blinded_df.attrs.get("blinding_audit_stats", {}))
    key_payload = {
        "created_utc": created,
        "method": "random_permutation",
        "generator": "random.Random(seed).shuffle over a fixed code-label pool",
        "seed": seed,
        "n_animals": len(codes),
        "codes": {str(t): codes[t] for t in sorted(codes)},
        "code_pool_size": len(_code_pool()),
        "blinded_manifest_sha256": _sha256(blinded_df.to_csv(index=False)),
        "audit_findings": findings,
        "audit_statistics": audit_stats,
        "custodian_note": (
            "This file IS the key. It must not be copied into the repository, into "
            "Dropbox, or onto shared storage, and it must not be given to whoever "
            "draws the regions, trains the classifiers or annotates validation data."),
    }
    _atomic_write(key_json, json.dumps(key_payload, indent=2, default=str) + "\n", 0o600)

    # ---- blinded output ---------------------------------------------------------
    _atomic_write(blinded_csv, blinded_df.to_csv(index=False), 0o644)
    report = {
        "created_utc": created,
        "n_rows": int(len(blinded_df)),
        "n_animals": int(blinded_df[_column(blinded_df, "code")].nunique()),
        "columns": list(map(str, blinded_df.columns)),
        "dropped_columns": split_report.get("dropped_columns", []),
        "notes": split_report.get("notes", []),
        "sorted_by": split_report.get("sorted_by", []),
        "audit_findings": findings,
        "audit_statistics": audit_stats,
        "n_leaks": len(hard),
        "n_risks": len([f for f in soft if f.startswith("RISK:")]),
        "seed": "NOT RECORDED HERE -- it lives with the key, in the custodian directory",
    }
    _atomic_write(report_json, json.dumps(report, indent=2, default=str) + "\n", 0o644)

    return {
        "custodian_dir": str(custodian),
        "work_dir": str(work),
        "key_csv": str(key_csv),
        "key_json": str(key_json),
        "provenance_csv": str(provenance_csv),
        "blinded_csv": str(blinded_csv),
        "report_json": str(report_json),
        "created_utc": created,
        "n_animals": len(codes),
        "n_private_rows": int(len(private_df)),
        "n_blinded_rows": int(len(blinded_df)),
        "blinded_columns": list(map(str, blinded_df.columns)),
        "dropped_columns": split_report.get("dropped_columns", []),
        "notes": split_report.get("notes", []),
        "audit": findings,
        "audit_statistics": audit_stats,
        "n_leaks": len(hard),
        "n_risks": len([f for f in soft if f.startswith("RISK:")]),
        "blinded_manifest_sha256": key_payload["blinded_manifest_sha256"],
    }
