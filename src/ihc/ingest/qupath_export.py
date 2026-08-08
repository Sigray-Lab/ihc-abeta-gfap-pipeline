"""Prepare a **blinded** QuPath project for manual region delineation.

Why this module exists
----------------------
Region delineation (spec §9) is done by hand in the QuPath GUI by the person who
also did the staining and the imaging. He can infer treatment group from a tube
number, so the project he opens must not contain one -- not in the file name, not
in the displayed image name, not in the slide label image, not in the internal
series name that QuPath shows in its metadata panel.

This module turns a *blinded manifest* into a JSON "project spec" that
``qupath/scripts/import_blinded_project.groovy`` reads to build that project. All
the work that can be done in Python -- resolving codes back to files, finding the
right series, checking the payload is present, auditing for identifier leaks --
is done here, where it is testable without a Java runtime.

Which sections go in, and why
-----------------------------
Every section whose staining condition is **resolved** goes in: the positives *and*
the negative controls.

*Positives* are obvious: the anatomical ROI is the denominator of every reported
percent-area number, so hippocampal formation and isocortex have to be outlined on
each positive section.

*Negative controls are not optional either.* ``config/config.yaml`` sets
``negative_control.assessed_per: [section, animal, channel, region]`` and the spec
(§7) requires reporting "classifier-positive area and background distribution in
every channel of every negative-control section". Both are **per region**. There is
no way to compute "Abeta-positive percent area inside hippocampus on a no-primary
section" without a hippocampus outline on that section, so a project containing only
positives silently makes the negative-control gate unevaluable. Delineating them
later, in a second pass, would also mean the ROIs on positives and negatives were
drawn at different times by a differently-calibrated hand, which is exactly the kind
of differential measurement error the spec spends §2 trying to avoid.

What is left out:

* Sections belonging to a slide flagged ``needs_confirmation``. Tube 37 is the live
  case: its ``positive_box`` column and its bench annotation name opposite boxes, so
  neither of its boxes has a resolved condition. ADR-0003 says such a row is not
  usable until it is confirmed at the bench; it is reported, never guessed. It is
  reported in ``spec["excluded"]``.
* Sections whose animal has no payload folder on disk. 23 of the 31 animals are
  index-only -- the metadata is known and the pixels have not been transferred yet.
  That is a normal state, not an error (spec §5), so they land in ``spec["skipped"]``
  with reason ``payload_absent`` and the run continues.
* The slide label series (``stack1``) and the slide overview, always. The label shows
  the tube ID both as printed text and as a DataMatrix barcode. Only series whose
  internal name ends in ``_01`` .. ``_04`` are ever referenced, which makes including
  the label a structural impossibility rather than something the caller must remember.

Leak vectors this module is written against
-------------------------------------------
All of these are present in this dataset and all of them are handled here:

===================================== =========================================
leak                                   handling
===================================== =========================================
``Image_29.vsi`` in the path           coded symlink ``K07.vsi`` -> the real file
``_Image_29_`` payload folder          coded symlink ``_K07_`` beside it
``stack1``, the slide label image      never referenced; only ``_0N`` series are
internal series name ``60_20x_...``    never written to the spec; the Groovy
                                       matches on the ``_0N`` suffix alone
displayed image name                   set explicitly to ``<code>_s<NN>``
acquisition timestamp (= group order)  never written to the spec
exposure (identifies tubes 51 and 60)  never written to the spec
project ordering = file order          entries sorted by code, not by tube
codes that preserve tube order         :func:`audit_codes` refuses them
===================================== =========================================

The one leak that cannot be closed from here is the *content* of the raw file: VSI
tag 2061 holds the tube ID and Bio-Formats will surface it to anyone who opens
QuPath's full metadata viewer. Nothing short of transcoding the pixels removes it.
It is documented for the delineator in ``docs/delineation_instructions.md`` as a
"do not go looking" instruction, and is stated as a residual limitation.

Public API
----------
:func:`build_project_spec` -- pure, read-only; returns the spec dict.
:func:`write_project_spec` -- writes the JSON and materialises the coded symlinks.
:func:`audit_codes` -- standalone check that a code mapping is not order-preserving.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

from ihc.ingest.verify import find_companion, read_vsi_index
from ihc.ingest.vsi_meta import VsiParseError, read_vsi_meta

__all__ = [
    "ProjectSpecError",
    "BlindingLeakError",
    "SPEC_SCHEMA",
    "build_project_spec",
    "write_project_spec",
    "audit_codes",
    "scan_for_tube_identifiers",
]

SPEC_SCHEMA = "ihc.ingest.qupath_export/1"

#: Pinned in ``env/tool_versions.yaml`` (D-15). Recorded in the spec so the Groovy
#: can warn when it is running under a different build.
DEFAULT_QUPATH_VERSION = "0.7.0"

#: Tube numbers that exist in this cohort. 59 was excluded before imaging.
COHORT_TUBE_IDS = frozenset(list(range(29, 59)) + [60])

#: How the scanner names things in ``raw_root``.
DEFAULT_VSI_NAME_PATTERN = "Image_{tube_id}.vsi"
DEFAULT_RESCAN_SUBDIR = "Rescan"

#: Only these series are ever opened. Anchored on the trailing ``_0N`` because the
#: series names are NOT uniform -- tube 60's read ``60_20x_DAPI, FITC, Cy3_01``.
SECTION_LABEL_RE = re.compile(r"^0[1-9]$")

#: Series names that must never be opened, checked again on the Groovy side.
FORBIDDEN_SERIES_NAMES = ("Label", "Overview", "Sample Mask", "FocusMap", "FocusPoints")

# --- identifier-leak patterns ---------------------------------------------- #

#: ``Image_29``, ``image 29``, ``Image-29``. Safe to run over a whole path.
_FILENAME_TUBE_RE = re.compile(r"(?i)image[ _\-]?\d{1,4}")

#: A bare cohort tube number standing on its own. The lookarounds are deliberately
#: alphanumeric rather than ``\b`` so that a coded ID like ``K42`` does NOT match
#: (there is no boundary between ``K`` and ``4``) while ``1007344 - 29``, the text
#: printed on the slide label, does.
_BARE_TUBE_RE = re.compile(r"(?<![0-9A-Za-z])(?:29|[3-5][0-9]|60)(?![0-9A-Za-z])")

#: A coded ID must look like this. Deliberately narrow: no spaces, no punctuation
#: that would need escaping in a file name, and a leading letter so it can never be
#: read as a number.
_CODE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]{1,15}$")

#: Column-name synonyms accepted in the blinded manifest, so this module does not
#: dictate the manifest's column names to the stage that builds it.
_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "coded_id": ("coded_id", "code", "blinded_id", "coded_animal_id"),
    "section_label": ("section_label", "section", "section_number"),
    # The physical section identity, stable across imaging sessions -- unlike
    # section_label, which each session renumbers from 01.
    "physical_section_label": ("physical_section_label", "physical_section"),
    "condition": ("condition", "stain_condition", "positive_negative", "pos_neg"),
    "scan": ("scan", "scan_source", "acquisition"),
    "needs_confirmation": ("needs_confirmation", "needs_confirm", "unresolved"),
    "series_index": ("series_index", "bioformats_series_index"),
}

_POSITIVE_WORDS = {"positive", "pos", "p", "primary"}
_NEGATIVE_WORDS = {"negative", "neg", "n", "control", "negative_control", "no_primary"}
#: ``manifest.py`` gives ``condition`` three states, not two. The third is what a
#: section gets when ``slides.csv`` has no row for it, when the row is flagged
#: ``needs_confirmation``, or when a rescan section could not be matched back to an
#: original. Spec §2: "If a section's condition cannot be resolved, exclude it and
#: report it" -- so it is an exclusion, never an exception and never a guess.
_UNRESOLVED_WORDS = {"unresolved", "unknown", "tbd", "pending", "na", "n_a"}


class ProjectSpecError(RuntimeError):
    """The project spec could not be built, or is not safe to write."""


class BlindingLeakError(ProjectSpecError):
    """A string destined for the blinded project carries an animal identifier.

    This is always a hard failure. A blinded artefact that leaks is worse than no
    blinded artefact at all, because the leak is invisible in the output and the
    resulting delineation looks exactly as trustworthy as a clean one.
    """


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #


def _is_blank(value: Any) -> bool:
    """True for ``None``, an empty/whitespace string, or a pandas ``NaN``."""
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    if isinstance(value, str):
        return value.strip() == "" or value.strip().lower() in {"nan", "none"}
    return False


def _rows_of(blinded_df: Any) -> list[dict[str, Any]]:
    """Normalise a DataFrame / list of mappings into a list of plain dicts."""
    if blinded_df is None:
        raise ProjectSpecError("blinded_df is None; expected a DataFrame or rows")
    to_dict = getattr(blinded_df, "to_dict", None)
    if callable(to_dict):  # pandas DataFrame
        try:
            return list(to_dict(orient="records"))
        except TypeError:  # a plain dict has .to_dict? no -- but be safe
            pass
    if isinstance(blinded_df, Mapping):
        raise ProjectSpecError(
            "blinded_df is a single mapping; expected a table of rows "
            "(a DataFrame, or an iterable of row mappings)"
        )
    if isinstance(blinded_df, Iterable):
        rows = []
        for i, row in enumerate(blinded_df):
            if not isinstance(row, Mapping):
                raise ProjectSpecError(f"row {i} is {type(row).__name__}, expected a mapping")
            rows.append(dict(row))
        return rows
    raise ProjectSpecError(f"cannot read rows from a {type(blinded_df).__name__}")


def _get(row: Mapping[str, Any], field: str) -> Any:
    """Read a logical field from a row, accepting the documented column synonyms."""
    for name in _COLUMN_ALIASES[field]:
        if name in row and not _is_blank(row[name]):
            return row[name]
    return None


def _normalise_condition(value: Any, context: str) -> str:
    """Return ``positive``, ``negative`` or ``unresolved``.

    ``unresolved`` is a legitimate third state, not an error: ``manifest.py`` assigns
    it to any section whose staining condition cannot be established from
    ``slides.csv``. The caller excludes those and reports them.
    """
    if _is_blank(value):
        return "unresolved"
    text = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if text in _POSITIVE_WORDS:
        return "positive"
    if text in _NEGATIVE_WORDS:
        return "negative"
    if text in _UNRESOLVED_WORDS:
        return "unresolved"
    raise ProjectSpecError(
        f"{context}: condition {value!r} is not one of positive / negative / "
        f"unresolved. Refusing to interpret it: a mislabelled condition means "
        f"negative controls are quantified as data and every downstream number "
        f"still looks entirely plausible (spec §2)."
    )


def _normalise_section_label(value: Any, context: str) -> str:
    text = str(value).strip()
    if text.isdigit() and len(text) == 1:
        text = f"0{text}"
    if not SECTION_LABEL_RE.match(text):
        raise ProjectSpecError(
            f"{context}: section label {value!r} is not of the form 01..09"
        )
    return text


# --------------------------------------------------------------------------- #
# identifier-leak auditing
# --------------------------------------------------------------------------- #


def scan_for_tube_identifiers(text: str, *, strict: bool = True) -> list[str]:
    """Return every animal-identifier-looking substring in *text*.

    Args:
        text: the string to inspect.
        strict: when True (the default) both a filename pattern (``Image_29``) and a
            bare cohort tube number (``29`` .. ``60`` standing alone) are reported.
            When False only the filename pattern is reported -- use this for strings
            that legitimately contain unrelated numbers, such as a timestamp.

    Returns:
        The matched substrings, in order of appearance. Empty means clean.
    """
    if not isinstance(text, str):
        return []
    hits = [m.group(0) for m in _FILENAME_TUBE_RE.finditer(text)]
    if strict:
        hits += [m.group(0) for m in _BARE_TUBE_RE.finditer(text)]
    return hits


def _assert_clean(text: str, what: str, *, strict: bool = True) -> None:
    hits = scan_for_tube_identifiers(text, strict=strict)
    if hits:
        raise BlindingLeakError(
            f"{what} contains what looks like an animal identifier "
            f"({', '.join(sorted(set(hits)))}): {text!r}\n"
            f"  This string would end up in the blinded QuPath project. Refusing to "
            f"build it. If the match is a false positive (an unrelated number in a "
            f"directory name), rename that directory rather than relaxing the check."
        )


def _assert_clean_path(text: str, what: str) -> None:
    """Identifier check for a path, which needs different rules from a plain name.

    The leak in a path is the *file name*: ``Image_29.vsi``, ``_Image_29_``. Parent
    directories are named by the custodian and their numbers mean nothing -- a run
    directory called ``pytest-33`` or ``batch-42`` is not an animal. Applying the
    bare-number rule to a whole path therefore produces false refusals, which is not
    a safe kind of noise: a check that cries wolf gets relaxed.

    So: the filename pattern (``Image_29``) runs over the whole path, because a
    directory named after a slide is a real leak; the bare-number pattern runs over
    the basename only.
    """
    _assert_clean(text, what, strict=False)
    _assert_clean(os.path.basename(text.rstrip("/")), f"{what} (file name)", strict=True)


def _redact(text: str) -> str:
    """Blank out anything identifier-shaped in free text bound for the blinded spec.

    Diagnostic messages are the sneakiest leak in the whole module. A warning raised
    deep inside the metadata reader quotes the internal series name, and tube 60's
    internal series names are ``60_20x_DAPI, FITC, Cy3_01`` -- so a message about a
    missing pixel size would carry the tube number into a "blinded" artefact. Paths in
    "file not found" messages do the same. Everything free-text is therefore redacted
    on the way in, and :func:`_assert_spec_clean` verifies the result.
    """
    text = _FILENAME_TUBE_RE.sub("<redacted>", str(text))
    return _BARE_TUBE_RE.sub("<redacted>", text)


def _walk_strings(node: Any, path: str = "") -> Iterable[tuple[str, str]]:
    """Yield ``(json_path, string)`` for every string value inside *node*."""
    if isinstance(node, str):
        yield path or "<root>", node
    elif isinstance(node, Mapping):
        for key, value in node.items():
            yield from _walk_strings(value, f"{path}.{key}" if path else str(key))
    elif isinstance(node, (list, tuple)):
        for i, value in enumerate(node):
            yield from _walk_strings(value, f"{path}[{i}]")


def _assert_spec_clean(spec: Mapping[str, Any], *, exempt_suffixes: tuple[str, ...] = ()) -> None:
    """Final gate: no string anywhere in the written spec may carry an identifier.

    This is the backstop. Every individual field is checked as it is built, but a
    field added later and forgotten here is exactly how a blinded artefact starts
    leaking, so the finished structure is walked once more before it is returned.

    ``link_plan`` is the only exempt subtree: it is the custodian's plan for building
    the coded symlinks and necessarily names the real files, and it is stripped before
    the spec is written to disk. Path-shaped fields use :func:`_assert_clean_path`.
    *exempt_suffixes* exempts named leaf fields, and is used for one thing only:
    ``vsi_path`` when the caller has explicitly and in writing accepted uncoded paths.
    """
    for key, value in spec.items():
        if key == "link_plan":
            continue
        for json_path, text in _walk_strings(value, key):
            if any(json_path.endswith(suffix) for suffix in exempt_suffixes):
                continue
            leaf = json_path.rsplit(".", 1)[-1]
            if leaf.endswith(("_path", "_dir", "_link")):
                _assert_clean_path(text, f"spec field {json_path}")
            else:
                _assert_clean(text, f"spec field {json_path}")


def audit_codes(codes: Mapping[Any, str]) -> list[str]:
    """Check that a tube-to-code mapping is a real blinding, and return warnings.

    Args:
        codes: ``{tube_id: coded_id}``.

    Returns:
        Non-fatal warnings, e.g. a code that happens to embed its own tube number.

    Raises:
        ProjectSpecError: the mapping is empty, has duplicate codes, has codes in a
            format that cannot be safely used as a file name, or -- the important
            one -- is **order-preserving**.

    Why order-preservation is fatal rather than a warning. Tube IDs run in
    contiguous treatment blocks (29-40 Rapamycin Diet, 41-48 Extra Control Diet,
    49-54 Control IP, 55-60 Rapamycin IP). Any code assignment that preserves the
    order of the tube IDs -- sequential numbering in file-iteration order, a hash
    that happens to sort the same way, ``code = tube + k`` -- reproduces the group
    structure exactly, so the "blinded" IDs still cluster by group and the blinding
    buys nothing (spec §2). Under a genuine random permutation of 31 animals the
    chance of the codes coming out monotonic is 1/31!, so this can be refused
    outright with no risk of a false alarm.
    """
    if not codes:
        raise ProjectSpecError("the code mapping is empty")

    normalised: dict[int, str] = {}
    for tube, code in codes.items():
        try:
            tube_int = int(str(tube).strip())
        except (TypeError, ValueError) as exc:
            raise ProjectSpecError(
                f"code mapping key {tube!r} is not a tube number. `codes` maps "
                f"tube_id -> coded_id; it looks like it may be the other way round."
            ) from exc
        code_text = str(code).strip()
        if not _CODE_RE.match(code_text):
            raise ProjectSpecError(
                f"coded ID {code!r} (tube {tube_int}) is not of the form "
                f"letter + alphanumerics, max 16 chars. Codes become file names and "
                f"displayed image names, so they must be boring."
            )
        normalised[tube_int] = code_text

    if len(set(normalised.values())) != len(normalised):
        seen: dict[str, int] = {}
        clashes = []
        for tube, code in sorted(normalised.items()):
            if code in seen:
                clashes.append(f"{code} (tubes {seen[code]} and {tube})")
            seen[code] = tube
        raise ProjectSpecError(f"duplicate coded IDs: {', '.join(clashes)}")

    warnings: list[str] = []
    outside = sorted(set(normalised) - COHORT_TUBE_IDS)
    if outside:
        warnings.append(
            f"the code mapping covers tube(s) {outside}, which are not in this cohort "
            f"(29-58 and 60; 59 was excluded before imaging). Check the mapping is the "
            f"right way round and belongs to this study."
        )
    for tube, code in sorted(normalised.items()):
        _assert_clean(code, f"coded ID for tube {tube}")
        if str(tube) in code:
            warnings.append(
                f"coded ID {code!r} contains the digits of its own tube number "
                f"({tube}). Not a leak on its own, but reseed if it worries you."
            )

    # Order preservation. Compare the numeric tails of the codes in tube order.
    ordered = [normalised[t] for t in sorted(normalised)]
    tails = []
    for code in ordered:
        digits = re.sub(r"\D", "", code)
        tails.append(int(digits) if digits else None)
    if len(ordered) >= 6 and all(t is not None for t in tails):
        if all(a < b for a, b in zip(tails, tails[1:])):
            raise ProjectSpecError(
                "the coded IDs increase monotonically with tube ID, so they preserve "
                "the tube order exactly. Tube IDs run in contiguous treatment blocks, "
                "so this reproduces the group structure and blinds nobody (spec §2). "
                "Draw the codes from a random permutation with a recorded seed."
            )
        if all(a > b for a, b in zip(tails, tails[1:])):
            raise ProjectSpecError(
                "the coded IDs decrease monotonically with tube ID. Reversing the "
                "order still preserves it. Use a random permutation (spec §2)."
            )
    if all(t is not None for t in tails):
        identical = [t for t, code_tail in zip(sorted(normalised), tails) if t == code_tail]
        if len(identical) > 1:
            warnings.append(
                f"{len(identical)} coded IDs carry the same number as their tube "
                f"(tubes {identical[:5]}{'...' if len(identical) > 5 else ''})"
            )
    return warnings


# --------------------------------------------------------------------------- #
# raw-file resolution
# --------------------------------------------------------------------------- #


def _stack_ids(companion: str) -> list[int]:
    """Return the ``stackNNNNN`` IDs inside a payload folder, ascending.

    Kept local rather than imported so this module owns its own contract: the
    ordering here is load-bearing (it is what the Bio-Formats series index is
    derived from) and must not drift if the verifier's private helper changes.
    """
    found: list[int] = []
    try:
        entries = os.listdir(companion)
    except OSError:
        return []
    for entry in entries:
        if entry.startswith("stack") and entry[5:].isdigit():
            if os.path.isdir(os.path.join(companion, entry)):
                found.append(int(entry[5:]))
    return sorted(found)


def _series_index_map(vsi_path: str, companion: str) -> tuple[dict[str, int], list[str]]:
    """Map ``section_label -> Bio-Formats series index``, plus any warnings.

    **This is a hint, not a guarantee, and the Groovy script re-derives it.**

    The ``.vsi`` index lists 15 records for a four-section slide (Label, Overview,
    Sample Mask, then each tissue series followed by its FocusMap and FocusPoints),
    but the payload holds only six stacks and Bio-Formats reports six series. So the
    position in the ``.vsi`` index is *not* the series index. What is derived here is
    the position of the section's stack folder in the ascending stack list --
    ``stack1`` (label) = 0, ``stack10000`` (overview) = 1, ``stack10002`` = 2, and so
    on -- which is the order Bio-Formats walks the payload in.

    Whether Bio-Formats 8.5.0 enumerates VSI series in exactly that order could not
    be checked on this machine (no Java, no QuPath). The import script therefore
    treats the index as advisory: it resolves the series by matching the internal
    name's trailing ``_0N``, cross-checks the pixel dimensions, and refuses rather
    than opening a series it cannot confirm.
    """
    warnings: list[str] = []
    stack_ids = _stack_ids(companion)
    if not stack_ids:
        return {}, [f"{os.path.basename(companion)}: no stackNNNNN folders"]

    try:
        index = read_vsi_index(vsi_path)
    except (OSError, ValueError) as exc:  # pragma: no cover - unreadable index
        return {}, [f"{os.path.basename(vsi_path)}: cannot read the index: {exc}"]

    tissue = [s for s in index if s.get("kind") == "tissue"]
    # Same rule the verifier uses: tissue series in index order bind to the tissue
    # stack folders in ascending stack order, which is acquisition order for both.
    tissue_stacks = [i for i in stack_ids if i not in (1, 10000)]
    if len(tissue) != len(tissue_stacks):
        warnings.append(
            f"{os.path.basename(vsi_path)}: the index names {len(tissue)} tissue "
            f"series but the payload holds {len(tissue_stacks)} tissue stack(s); "
            f"series indices left unset and resolved by name at import time"
        )
        return {}, warnings

    mapping: dict[str, int] = {}
    for series, stack_id in zip(tissue, tissue_stacks):
        label = series.get("section")
        if label:
            mapping[str(label)] = stack_ids.index(stack_id)
    return mapping, warnings


def _resolve_vsi_path(raw_root: Path, tube_id: int, scan: str, *,
                      vsi_name_pattern: str, rescan_subdir: str) -> Path:
    directory = raw_root if scan == "original" else raw_root / rescan_subdir
    return directory / vsi_name_pattern.format(tube_id=tube_id)


def _alias_stem(code: str, scan: str) -> str:
    """Coded file-name stem. A second scan of the same animal gets a ``b`` suffix."""
    return code if scan == "original" else f"{code}b"


# --------------------------------------------------------------------------- #
# channels
# --------------------------------------------------------------------------- #


def _channel_plan(channels_cfg: Mapping[str, Any] | None) -> dict[str, Any]:
    """Decide which channels exist and which one stays visible.

    Only the nuclear counterstain is visible. The marker channels are hidden because
    anatomical boundaries must be drawn without sight of plaque or GFAP signal
    (spec §9; ``regions.registration.hide_marker_channels_while_delineating: true``).
    This is a scientific requirement: a boundary nudged towards or away from visible
    burden biases the denominator of every number computed inside it, and it does so
    even though the person drawing it is blinded to treatment group.
    """
    if channels_cfg is None:
        from ihc.util.config import load_channels  # local import: keeps this testable

        channels_cfg = load_channels()

    order = list(channels_cfg.get("expected_order") or [])
    entries = list(channels_cfg.get("channels") or [])
    if not order and entries:
        order = [str(c.get("name")) for c in entries]
    if not order:
        raise ProjectSpecError("config/channels.yaml declares no channel order")

    targets = {str(c.get("name")): str(c.get("target", "")) for c in entries}
    visible = [name for name in order if targets.get(name) == "nuclei"]
    if not visible:
        visible = [name for name in order if name.upper() == "DAPI"]
    if len(visible) != 1:
        raise ProjectSpecError(
            f"expected exactly one nuclear-counterstain channel to leave visible, "
            f"found {visible!r} in config/channels.yaml"
        )
    hidden = [name for name in order if name not in visible]
    return {
        "names": order,
        "visible": visible,
        "hidden": hidden,
        "reason_hidden": (
            "Anatomical boundaries must be drawn without seeing marker signal, or "
            "visible plaque/GFAP burden influences where the boundary goes "
            "(CLAUDE_v1.2 §9)."
        ),
    }


# --------------------------------------------------------------------------- #
# the spec
# --------------------------------------------------------------------------- #


def build_project_spec(
    blinded_df: Any,
    codes: Mapping[Any, str],
    raw_root: str | os.PathLike[str],
    out_dir: str | os.PathLike[str],
    *,
    project_name: str = "delineation",
    include_conditions: tuple[str, ...] = ("positive", "negative"),
    alias_images: bool = True,
    acknowledge_uncoded_paths: str | None = None,
    channels_cfg: Mapping[str, Any] | None = None,
    qupath_version: str = DEFAULT_QUPATH_VERSION,
    vsi_name_pattern: str = DEFAULT_VSI_NAME_PATTERN,
    rescan_subdir: str = DEFAULT_RESCAN_SUBDIR,
) -> dict:
    """Build the JSON-serialisable spec describing a blinded QuPath project.

    Read-only: nothing is written and no symlink is created here. Pass the result to
    :func:`write_project_spec`, which materialises it.

    Args:
        blinded_df: the blinded analysis manifest -- a pandas DataFrame, or any
            iterable of row mappings. One row per section. Required fields (column
            synonyms in ``_COLUMN_ALIASES`` are accepted): ``coded_id``,
            ``section_label``, ``condition``. Optional: ``scan``
            (``original`` | ``rescan``, default ``original``), ``needs_confirmation``,
            ``series_index`` (cross-checked against the one derived here).
            The manifest carries **no** tube ID -- that is the point of it -- so the
            join back to the raw files goes through *codes*. Note that
            ``blinding.NEVER_BLINDED`` blocks any column whose name contains "scan",
            so a blinded manifest will not carry one and every row is read as the
            original acquisition. Selecting the rescan of tubes 51 and 60 (which
            ``config.yaml`` prefers) is therefore the custodian's to add, by setting
            ``scan="rescan"`` on those rows before calling this.
        codes: ``{tube_id: coded_id}``, from the custodian's blinding key. Audited by
            :func:`audit_codes` before anything else happens.
        raw_root: the read-only raw-data root holding ``Image_NN.vsi`` and the
            ``_Image_NN_`` payload folders, with rescans in ``<raw_root>/Rescan``.
        out_dir: where the project will be built. The coded symlinks go in
            ``<out_dir>/images`` and the QuPath project in ``<out_dir>/qupath``.
        project_name: display name of the project. Must be identifier-free.
        include_conditions: which staining conditions to include. The default is
            **both**, and the module docstring explains why the negative controls are
            needed: the negative-control QC gate is assessed per region, so a
            negative section without an ROI cannot be evaluated at all.
        alias_images: when True (default), the spec references coded symlinks so that
            no path in the blinded project contains a tube ID. Set False only with
            *acknowledge_uncoded_paths*, and expect the Groovy script to refuse.
        acknowledge_uncoded_paths: a free-text acknowledgement (who accepted the risk
            and why) required when ``alias_images=False``. It is recorded in the spec
            and is what allows the import script to proceed with raw paths.
        channels_cfg: parsed ``config/channels.yaml``; loaded from the repo if None.
        qupath_version: the pinned version (D-15), recorded for the import script.
        vsi_name_pattern: how a tube number becomes a file name in ``raw_root``.
        rescan_subdir: subfolder of ``raw_root`` holding the re-acquired slides.

    Returns:
        A JSON-serialisable dict. The interesting keys:

        * ``images`` -- one entry per section to import, **sorted by coded name** so
          that project order carries no acquisition-order information. Each entry has
          ``vsi_path`` (absolute), ``series_index`` (advisory), ``series_match_suffix``
          (authoritative), ``image_name`` (``<code>_s<NN>``), ``channel_names``,
          ``channel_visible`` and the expected pixel dimensions.
        * ``skipped`` -- sections that could not be included but are not errors,
          each with a ``reason`` (``payload_absent`` is the common one).
        * ``excluded`` -- sections deliberately left out, e.g. ``needs_confirmation``.
        * ``warnings`` -- everything a human should read before handing the project on.

    Raises:
        ProjectSpecError: the manifest or the code mapping is unusable.
        BlindingLeakError: a string bound for the project carries an identifier.

    Example:
        >>> spec = build_project_spec(df, {29: "K07"}, raw, out)   # doctest: +SKIP
        >>> spec["images"][0]["image_name"]                        # doctest: +SKIP
        'K07_s01'
    """
    warnings: list[str] = list(audit_codes(codes))

    if alias_images and acknowledge_uncoded_paths:
        raise ProjectSpecError(
            "acknowledge_uncoded_paths was given but alias_images is True. "
            "The acknowledgement only applies when aliasing is switched off."
        )
    if not alias_images and not (acknowledge_uncoded_paths or "").strip():
        raise ProjectSpecError(
            "alias_images=False makes the raw .vsi paths -- which contain the tube "
            "number -- visible inside the blinded QuPath project. If that is a "
            "considered decision, record it: pass acknowledge_uncoded_paths='<who, "
            "why>'. The import script will still report every path it finds."
        )

    raw_root = Path(raw_root).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    images_dir = out_dir / "images"

    _assert_clean(project_name, "project_name")
    _assert_clean_path(str(out_dir), "out_dir")

    code_to_tube: dict[str, int] = {}
    for tube, code in codes.items():
        code_to_tube[str(code).strip()] = int(str(tube).strip())

    channel_plan = _channel_plan(channels_cfg)
    for name in channel_plan["names"]:
        _assert_clean(str(name), "channel name")

    rows = _rows_of(blinded_df)
    if not rows:
        raise ProjectSpecError("the blinded manifest has no rows")

    images: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    link_plan: list[dict[str, Any]] = []

    # Cache per (tube, scan): resolving a .vsi costs a full index parse.
    slide_cache: dict[tuple[int, str], dict[str, Any]] = {}

    for i, row in enumerate(rows):
        context = f"manifest row {i}"

        code = _get(row, "coded_id")
        if code is None:
            raise ProjectSpecError(
                f"{context}: no coded ID. Expected one of "
                f"{_COLUMN_ALIASES['coded_id']}."
            )
        code = str(code).strip()
        if code not in code_to_tube:
            raise ProjectSpecError(
                f"{context}: coded ID {code!r} is not in the code mapping, so it "
                f"cannot be resolved to a file. The manifest and the blinding key "
                f"disagree; do not guess."
            )
        tube_id = code_to_tube[code]

        section_label = _normalise_section_label(_get(row, "section_label"), context)
        scan_raw = _get(row, "scan")
        scan = str(scan_raw).strip().lower() if scan_raw is not None else "original"
        if scan in {"", "original", "first", "main"}:
            scan = "original"
        elif scan in {"rescan", "redo", "second"}:
            scan = "rescan"
        else:
            raise ProjectSpecError(f"{context}: unknown scan {scan_raw!r}")

        entry_stub = {"code": code, "section_label": section_label, "scan": scan,
                      "physical_section_label": _get(row, "physical_section_label")}

        needs_confirmation = _get(row, "needs_confirmation")
        if not _is_blank(needs_confirmation):
            excluded.append({
                **entry_stub,
                "reason": "needs_confirmation",
                "detail": (
                    "the slide record is flagged as unresolved, so this section's "
                    "staining condition is not known (ADR-0003). Excluded and "
                    "reported, never guessed."
                ),
            })
            continue

        condition = _normalise_condition(_get(row, "condition"), context)
        if condition == "unresolved":
            excluded.append({
                **entry_stub,
                "reason": "condition_unresolved",
                "detail": (
                    "the staining condition could not be established from the slide "
                    "record, so this section is neither a datum nor a control. "
                    "Excluded and reported (spec §2)."
                ),
            })
            continue
        if condition not in include_conditions:
            excluded.append({**entry_stub, "reason": f"condition_{condition}_not_requested"})
            continue

        # `<code>_s<NN>`, where NN is the PHYSICAL section identity, not the scanner's
        # label for it. The two differ: a second imaging session numbers its sections
        # from 01 regardless of where they sit on the slide, so the scanner label
        # collides across scans while denoting different pieces of tissue. Keying the
        # display name on it merged sections that are not the same section, and the
        # duplicate guard below then excluded both. Physical identity is stable across
        # scans and is what the manifest already resolves.
        alias_stem = _alias_stem(code, scan)
        physical = str(entry_stub.get("physical_section_label") or section_label)
        image_name = f"{alias_stem}_s{physical}"
        _assert_clean(image_name, f"{context}: displayed image name")

        key = (tube_id, scan)
        if key not in slide_cache:
            slide_cache[key] = _inspect_slide(
                raw_root, tube_id, scan,
                vsi_name_pattern=vsi_name_pattern, rescan_subdir=rescan_subdir,
            )
        slide = slide_cache[key]

        if slide["skip_reason"] is not None:
            skipped.append({**entry_stub, "reason": slide["skip_reason"],
                            "detail": _redact(slide["skip_detail"])})
            continue

        series = slide["sections"].get(section_label)
        if series is None:
            skipped.append({
                **entry_stub,
                "reason": "section_not_in_file",
                "detail": (
                    f"the slide holds sections "
                    f"{sorted(slide['sections'])} and not {section_label}"
                ),
            })
            continue

        if alias_images:
            vsi_path = images_dir / f"{alias_stem}.vsi"
            payload_alias = images_dir / f"_{alias_stem}_"
        else:
            vsi_path = Path(slide["vsi_path"])
            payload_alias = None

        series_index = slide["series_index"].get(section_label)
        declared = _get(row, "series_index")
        if declared is not None and series_index is not None:
            if int(declared) != int(series_index):
                warnings.append(
                    f"{image_name}: the manifest declares series_index "
                    f"{int(declared)} but the payload layout implies "
                    f"{series_index}. The import script resolves by series name, so "
                    f"neither number is trusted blindly -- but check why they differ."
                )

        entry = {
            "image_name": image_name,
            "code": code,
            "section_label": section_label,
            "scan": scan,
            "condition": condition,
            "vsi_path": str(vsi_path),
            "series_index": series_index,
            "series_match_suffix": f"_{section_label}",
            "n_series_expected": slide["n_stacks"],
            "channel_names": list(channel_plan["names"]),
            "channel_visible": [n in channel_plan["visible"] for n in channel_plan["names"]],
            "width_px": series["width_px"],
            "height_px": series["height_px"],
            "pixel_size_um": series["pixel_size_um"],
        }
        if alias_images:
            entry["alias"] = {
                "stem": alias_stem,
                "vsi_link": str(vsi_path),
                "payload_link": str(payload_alias),
            }
            plan = {
                "stem": alias_stem,
                "vsi_link": str(vsi_path),
                "payload_link": str(payload_alias),
                "target_vsi": slide["vsi_path"],
                "target_payload": slide["companion"],
            }
            if plan not in link_plan:
                link_plan.append(plan)

        for field in ("image_name", "code", "series_match_suffix"):
            _assert_clean(str(entry[field]), f"{image_name}: {field}")
        if alias_images:
            # With aliasing off the path deliberately names the real file; that is
            # what the written acknowledgement buys, and the import script reports
            # every such path in full before it proceeds.
            _assert_clean_path(str(entry["vsi_path"]), f"{image_name}: vsi_path")

        images.append(entry)
        for message in slide["warnings"]:
            message = f"{code}: {_redact(message)}"
            if message not in warnings:
                warnings.append(message)

    # --- collisions ---------------------------------------------------------
    # Two rows resolving to one image is not a duplicate to be de-duplicated; it is
    # an ambiguity. It has one expected cause and it is structural: tubes 51 and 60
    # were scanned twice, and `blinding.NEVER_BLINDED` drops every column whose name
    # contains "scan", so both scans arrive in the blinded manifest looking identical.
    #
    # Resolving it by preferring the rescan would be a guess of exactly the kind
    # manifest.py warns about: the rescan holds only two sections, and its `_01` is
    # not necessarily the original's `_01` -- the physical sections are matched by
    # stage X, not by label. Picking wrong binds a different coronal level, and
    # nothing downstream would notice. So the colliding entries are excluded together
    # and reported, and the custodian disambiguates by setting `scan` explicitly.
    by_name: dict[str, list[dict[str, Any]]] = {}
    for entry in images:
        by_name.setdefault(entry["image_name"], []).append(entry)
    images = []
    for name, group in sorted(by_name.items()):
        if len(group) == 1:
            images.append(group[0])
            continue
        conditions = sorted({e["condition"] for e in group})
        contradiction = (
            f" The rows also DISAGREE about the staining condition ({', '.join(conditions)}), "
            f"so they are not two copies of one section -- they denote different "
            f"physical sections that happen to share a label."
            if len(conditions) > 1 else ""
        )
        excluded.append({
            "code": group[0]["code"],
            "section_label": group[0]["section_label"],
            "scan": group[0]["scan"],
            "reason": "ambiguous_duplicate_rows",
            "conditions_seen": conditions,
            "detail": (
                f"{len(group)} manifest rows resolve to the image {name}, so it is not "
                f"determined which physical section should be delineated.{contradiction} "
                f"The expected cause is a slide that was scanned twice: the blinded "
                f"manifest drops every column whose name contains 'scan', so the two "
                f"acquisitions arrive indistinguishable even though the rescan covers "
                f"only one PAP-pen box and its section labels are its own acquisition "
                f"order, not the original's. Fix upstream by letting a scan "
                f"discriminator through to the blinded manifest, or set scan='rescan' "
                f"on the re-acquired rows before calling. Not guessed here: choosing "
                f"wrong would delineate the positive box while the record says "
                f"negative control."
            ),
        })
        warnings.append(
            f"{name}: EXCLUDED -- {len(group)} manifest rows resolve to it"
            + (f" with conflicting conditions {conditions}" if len(conditions) > 1 else "")
            + ". See spec['excluded']."
        )

    # Sorting by coded name is not cosmetic. Left in manifest order, the project
    # browser would list the animals in tube order, which is treatment-group order.
    images.sort(key=lambda e: (e["code"], e["scan"], e["section_label"]))

    # Drop links for slides that lost every section to a collision, so the images
    # directory mirrors exactly what is imported and nothing dangles.
    used_stems = {e["alias"]["stem"] for e in images if e.get("alias")}
    link_plan = [p for p in link_plan if p["stem"] in used_stems]

    n_positive = sum(1 for e in images if e["condition"] == "positive")
    spec = {
        "schema": SPEC_SCHEMA,
        "project_name": project_name,
        "out_dir": str(out_dir),
        "images_dir": str(images_dir) if alias_images else None,
        "project_dir": str(out_dir / "qupath"),
        "image_type": "FLUORESCENCE",
        "channels": channel_plan,
        "blinding": {
            "paths_are_coded": bool(alias_images),
            "uncoded_paths_acknowledged_by": acknowledge_uncoded_paths,
            "forbidden_series_names": list(FORBIDDEN_SERIES_NAMES),
            "note": (
                "Only series whose internal name ends in _0N are ever opened, so the "
                "slide label image (which shows the tube number as text and as a "
                "DataMatrix barcode) cannot be imported by construction."
            ),
        },
        "qupath": {
            "expected_version": qupath_version,
            "set_image_name_from_spec": True,
            "resolve_series_by_name_suffix": True,
        },
        "counts": {
            "images": len(images),
            "positive": n_positive,
            "negative": len(images) - n_positive,
            "animals": len({e["code"] for e in images}),
            "skipped": len(skipped),
            "excluded": len(excluded),
        },
        "images": images,
        "skipped": skipped,
        "excluded": excluded,
        "warnings": warnings,
        # Custodian-only: the real files each coded symlink points at. Stripped by
        # write_project_spec before anything is written to disk, and exempt from the
        # identifier scan for exactly that reason.
        "link_plan": link_plan,
        # DATE only, deliberately. A full ISO timestamp is a poor thing to put in a
        # blinded artefact twice over: `09:34:12` reads as a bare tube number to any
        # identifier scanner (ours and other people's), and acquisition-adjacent
        # timing is the leak vector this cohort is most exposed to, since scanning
        # ran in ascending tube order and therefore ascending group order. The exact
        # generation time is in the custodian's own records, where it belongs.
        "provenance": {
            "generator": "ihc.ingest.qupath_export.build_project_spec",
            "generated_utc_date": _dt.datetime.now(_dt.timezone.utc).date().isoformat(),
            "raw_root_is_recorded_here": False,
        },
    }
    _assert_spec_clean(spec, exempt_suffixes=() if alias_images else (".vsi_path",))
    return spec


def _inspect_slide(
    raw_root: Path,
    tube_id: int,
    scan: str,
    *,
    vsi_name_pattern: str,
    rescan_subdir: str,
) -> dict[str, Any]:
    """Resolve one animal's ``.vsi`` + payload and describe its tissue series."""
    result: dict[str, Any] = {
        "vsi_path": None,
        "companion": None,
        "sections": {},
        "series_index": {},
        "n_stacks": None,
        "warnings": [],
        "skip_reason": None,
        "skip_detail": None,
    }
    vsi_path = _resolve_vsi_path(
        raw_root, tube_id, scan,
        vsi_name_pattern=vsi_name_pattern, rescan_subdir=rescan_subdir,
    )
    if not vsi_path.is_file():
        result["skip_reason"] = "vsi_absent"
        result["skip_detail"] = f"no index file at {vsi_path}"
        return result
    result["vsi_path"] = str(vsi_path)

    companion = find_companion(str(vsi_path))
    if companion is None:
        # 23 of 31 animals are index-only. Metadata known, pixels not transferred.
        # A normal state, not an error (spec §5).
        result["skip_reason"] = "payload_absent"
        result["skip_detail"] = (
            "the .vsi index is present but its companion payload folder is not, so "
            "there are no pixels to open. Transfer at the parent-folder level, then "
            "rebuild the project."
        )
        return result
    result["companion"] = str(Path(companion).resolve())

    try:
        meta = read_vsi_meta(str(vsi_path))
    except VsiParseError as exc:
        result["skip_reason"] = "vsi_unreadable"
        result["skip_detail"] = str(exc)
        return result

    for message in meta.warnings:
        result["warnings"].append(f"tube-index metadata: {message}")

    for series in meta.series:
        result["sections"][series.section_label] = {
            "width_px": series.width_px,
            "height_px": series.height_px,
            "pixel_size_um": series.pixel_size_um,
        }

    index_map, index_warnings = _series_index_map(str(vsi_path), companion)
    result["series_index"] = index_map
    result["warnings"].extend(index_warnings)
    result["n_stacks"] = len(_stack_ids(companion)) or None
    return result


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #


def write_project_spec(spec: Mapping[str, Any], out_dir: str | os.PathLike[str]) -> Path:
    """Write the spec JSON and materialise the coded image links.

    Args:
        spec: the dict from :func:`build_project_spec`.
        out_dir: must match ``spec["out_dir"]``. Passed explicitly so a spec cannot
            be written somewhere it does not describe.

    Returns:
        Path to the written ``project_spec.json``.

    Raises:
        ProjectSpecError: *out_dir* disagrees with the spec, a link target is
            missing, or a link name is already taken by something else.

    The coded links are the mechanism that keeps the tube number out of the project.
    For each animal one symlink ``<code>.vsi`` points at ``Image_NN.vsi`` and one
    symlink ``_<code>_`` points at ``_Image_NN_``. Bio-Formats finds a VSI's pixel
    payload by looking for a sibling folder named ``_<file stem>_``, so the pair has
    to be renamed together for the alias to resolve.

    **Untested assumption, stated plainly:** that Bio-Formats' companion-folder lookup
    is satisfied by this pair of symlinks could not be verified here -- QuPath 0.7.0
    is not installed on this machine and neither is a Java runtime. Everything Python
    can check *is* checked (the links exist, resolve, and point at a folder holding
    the expected stacks). If the assumption turns out to be wrong, the import fails
    loudly on the first image rather than producing something subtly wrong, and the
    documented fallback is ``alias_images=False`` plus an explicit acknowledgement.
    """
    out_dir = Path(out_dir).expanduser().resolve()
    declared = spec.get("out_dir")
    if declared and Path(declared) != out_dir:
        raise ProjectSpecError(
            f"this spec describes {declared}, not {out_dir}. Rebuild it for the "
            f"directory you mean to write, rather than moving it after the fact."
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    Path(spec["project_dir"]).mkdir(parents=True, exist_ok=True)

    if spec.get("blinding", {}).get("paths_are_coded"):
        images_dir = Path(spec["images_dir"])
        images_dir.mkdir(parents=True, exist_ok=True)
        aliased = {e["alias"]["stem"] for e in spec["images"] if e.get("alias")}
        planned = {p["stem"] for p in spec.get("link_plan", [])}
        if aliased - planned:
            raise ProjectSpecError(
                f"the spec says paths are coded but there is no link plan for "
                f"{sorted(aliased - planned)}. Rebuild the spec; do not hand-edit it."
            )
        for plan in spec.get("link_plan", []):
            _link(Path(plan["vsi_link"]), Path(plan["target_vsi"]), "index file")
            _link(Path(plan["payload_link"]), Path(plan["target_payload"]), "payload folder")

    # The link plan names the real files, so it never reaches disk inside the
    # project directory. It exists only long enough to create the symlinks.
    on_disk = {k: v for k, v in spec.items() if k != "link_plan"}
    on_disk["blinding"] = dict(on_disk.get("blinding", {}))
    on_disk["blinding"]["link_plan_redacted"] = "link_plan" in spec
    _assert_spec_clean(
        on_disk,
        exempt_suffixes=() if on_disk["blinding"].get("paths_are_coded") else (".vsi_path",),
    )

    spec_path = out_dir / "project_spec.json"
    tmp_path = spec_path.with_suffix(".json.tmp")
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(on_disk, handle, indent=2, sort_keys=False)
        handle.write("\n")
    os.replace(tmp_path, spec_path)  # atomic: never a half-written spec
    return spec_path


def _link(link: Path, target: Path, what: str) -> None:
    """Create ``link -> target``, idempotently, refusing to overwrite anything else."""
    if not target.exists():
        raise ProjectSpecError(f"cannot link the {what}: {target} does not exist")
    if link.is_symlink():
        current = Path(os.readlink(link))
        if current == target:
            return
        raise ProjectSpecError(
            f"{link} already points at {current}, not {target}. A coded name has "
            f"been reused for a different animal -- refusing to repoint it, because "
            f"annotations already drawn under that name would silently move to "
            f"another mouse. Build the project in a fresh directory."
        )
    if link.exists():
        raise ProjectSpecError(f"{link} exists and is not a symlink; refusing to replace it")
    link.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(target, link)
