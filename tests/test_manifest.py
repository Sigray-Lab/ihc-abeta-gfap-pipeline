"""Tests for `ihc.ingest.manifest` -- the row-per-section analysis manifest.

WHAT CLASS OF BUG THIS FILE PROTECTS AGAINST
============================================
**Quantifying negative-control sections as data.**

Every slide carries three or four sections from one animal in two PAP-pen boxes.
The two sections in *one box* share a staining condition: one box got primary
antibody, the other got DAPI + secondary only.  Which box is the positive one
varies per slide and is a bench fact recorded in ``config/slides.csv``.

If the manifest gets that assignment wrong, the pipeline measures blank tissue
as if it were signal.  Nothing downstream can detect it: percent-area numbers
stay in a plausible range, the QC plots look fine, the classifier is unchanged,
and the effect is a dilution towards the null in whichever direction the error
happens to fall.  There is no image, no plot and no summary statistic that shows
it.  Only these tests do.

The specific errors guarded here, each of which has a live counter-example in
this cohort:

* ``_01``/``_02`` and ``_03``/``_04`` assumed to pair.  **False.**  Section
  number is *acquisition order*, not slide position -- tube 49 runs
  ``02 03 | 01 04`` in stage order and tube 45 runs ``03 04 | 01 02``.
* two positives and two negatives assumed universal.  **False.**  Layouts
  2+2, 2+1, 1+2 and ``both`` (no negative control at all) all occur.
* four sections assumed per slide.  **False** for tubes 30, 33, 34, 42, 53, 54.
* a contradicted record silently guessed rather than excluded (tube 37).
* absent pixels treated as an error rather than a normal state -- 23 of 31
  animals are index-only right now, and the manifest must still describe them.

Ground truth used below was read out of the 31 real ``.vsi`` index files and
cross-checked against ``config/slides.csv``; it is hard-coded on purpose, so
that a change in ``assign_boxes`` cannot quietly move the expectation with it.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from _stage3_helpers import (  # noqa: E402 - tests dir is on sys.path via pytest
    ALL_TUBES,
    BOTH_BOXES_TUBES,
    NEEDS_CONFIRMATION_TUBES,
    PAYLOAD_TUBES,
    THREE_SECTION_TUBES,
    as_bool,
    boxes_by_section,
    col,
    conditions_by_section,
    find_column,
    read_slides_csv,
    require_cohort_index_files,
    rescan_mask,
    resolve_column,
    rows_for_tube,
    section_key,
)

pd = pytest.importorskip("pandas")


# ==========================================================================
# fixtures
# ==========================================================================
@pytest.fixture(scope="module")
def manifest_mod():
    """`ihc.ingest.manifest`, or a clean skip while it is still being written."""
    return pytest.importorskip(
        "ihc.ingest.manifest",
        reason="ihc.ingest.manifest is not written yet (stage 3, in progress)",
    )


@pytest.fixture(scope="module")
def manifest(manifest_mod, data_root):
    """The real manifest, built once from all 31 index files.

    Built with ``include_rescans=True`` because that is the production default
    and the tests must see whatever extra rows it adds; every per-animal
    assertion filters rescans out explicitly via `rows_for_tube`.
    """
    require_cohort_index_files(data_root)
    df = manifest_mod.build_manifest(raw_root=data_root, include_rescans=True)
    assert isinstance(df, pd.DataFrame), f"build_manifest returned {type(df)}"
    assert len(df) > 0, "build_manifest produced an empty manifest"
    return df


@pytest.fixture(scope="module")
def slides() -> dict:
    return read_slides_csv()


# ==========================================================================
# shape and coverage
# ==========================================================================
@pytest.mark.requires_data
def test_every_tube_in_slides_csv_appears_and_no_extras(manifest, slides):
    """The manifest describes the cohort exactly -- no drop-outs, no strays.

    An animal missing from the manifest vanishes from the analysis without a
    single error message; an extra one (a rescan mistaken for a new animal, or a
    stray file) turns up as a phantom replicate.
    """
    in_manifest = set(int(t) for t in col(manifest, "tube_id").unique())
    in_csv = set(slides)
    assert in_manifest == in_csv, (
        f"missing from manifest: {sorted(in_csv - in_manifest)}; "
        f"unexpected in manifest: {sorted(in_manifest - in_csv)}"
    )
    assert in_csv == set(ALL_TUBES), "slides.csv itself no longer matches the cohort"


@pytest.mark.requires_data
@pytest.mark.parametrize("tube", sorted(THREE_SECTION_TUBES))
def test_three_section_slides_produce_three_rows(manifest, tube):
    """Six animals have three sections. A hard-coded 4 invents a phantom row.

    A phantom row is worse than a missing one: it carries a condition, so it
    contributes a measurement, and area-weighted aggregation would silently
    give that animal a section that does not exist.
    """
    rows = rows_for_tube(manifest, tube)
    assert len(rows) == 3, (
        f"tube {tube} has three tissue series but the manifest holds "
        f"{len(rows)} rows: {sorted(col(rows, 'section_label').map(section_key))}"
    )


@pytest.mark.requires_data
def test_four_section_slides_produce_four_rows(manifest):
    """The complement: every other animal contributes exactly four rows."""
    wrong = {}
    for tube in ALL_TUBES:
        if tube in THREE_SECTION_TUBES:
            continue
        count = len(rows_for_tube(manifest, tube))
        if count != 4:
            wrong[tube] = count
    assert not wrong, f"tubes with an unexpected section count: {wrong}"


@pytest.mark.requires_data
def test_section_labels_are_strings_not_integers(manifest):
    """`"01"` must not become `1`.

    A section number that silently becomes an integer breaks every join against
    slides.csv and the QuPath project, and sorts ``10`` before ``2``.  The spec
    calls this out by name (section 6, step 2).
    """
    labels = col(manifest, "section_label").tolist()
    offenders = [v for v in labels if not isinstance(v, str)]
    assert not offenders, (
        f"{len(offenders)} section labels are not strings, e.g. {offenders[:5]!r}. "
        "Zero-padded strings are required."
    )
    assert all(section_key(v) in {"01", "02", "03", "04"} for v in labels)


# ==========================================================================
# condition assignment, per layout
# ==========================================================================
@pytest.mark.requires_data
def test_layout_2_plus_2_near_label_positive(manifest):
    """Tube 29: stage order ``01 02 | 03 04``, positive box = near_label."""
    assert conditions_by_section(manifest, 29) == {
        "01": "positive",
        "02": "positive",
        "03": "negative",
        "04": "negative",
    }


@pytest.mark.requires_data
def test_layout_2_plus_1_near_label_positive(manifest):
    """Tube 30: three sections, ``01 02 | 03``, positive box = near_label."""
    assert conditions_by_section(manifest, 30) == {
        "01": "positive",
        "02": "positive",
        "03": "negative",
    }


@pytest.mark.requires_data
def test_layout_2_plus_1_far_label_positive_gives_a_single_positive(manifest):
    """Tube 42: ``01 02 | 03`` with the FAR box positive -- one positive section.

    The live case for ``sections.flag_animals_with_single_positive_section``.
    Getting the box backwards here doubles the positive area and inverts which
    sections are the control.
    """
    assert conditions_by_section(manifest, 42) == {
        "01": "negative",
        "02": "negative",
        "03": "positive",
    }


@pytest.mark.requires_data
@pytest.mark.parametrize("tube", [33, 54])
def test_layout_1_plus_2_near_label_positive(manifest, tube):
    """Tubes 33 and 54: the 1+2 layout -- one positive, two negatives.

    The near box holds a single section here.  An implementation that assumes
    the near box always holds two would mark ``02`` positive as well.
    """
    assert conditions_by_section(manifest, tube) == {
        "01": "positive",
        "02": "negative",
        "03": "negative",
    }


# ==========================================================================
# the two non-monotonic slides -- the headline regressions
# ==========================================================================
@pytest.mark.requires_data
def test_tube_49_regression_far_box_is_sections_01_and_04(manifest):
    """TUBE 49 REGRESSION -- the slide that refutes the ``01+02 / 03+04`` rule.

    Sorted by stage X, tube 49 runs ``_02 _03`` -- 14.65 mm gap -- ``_01 _04``.
    The operator scanned the physically-third section first, so the section
    number is acquisition order, not slide position.  ``slides.csv`` says
    ``positive_box=far_label``, and the wet lab confirmed in writing on 2026-08-06 that
    the far/right box is the one that got primary antibody.

    Therefore: **01 and 04 are positive; 02 and 03 are the negative controls.**

    Any implementation that pairs on the section number gets exactly the
    opposite answer for three of the four sections, and every number this animal
    contributes is then measured off the wrong tissue.
    """
    assert conditions_by_section(manifest, 49) == {
        "01": "positive",
        "04": "positive",
        "02": "negative",
        "03": "negative",
    }


@pytest.mark.requires_data
def test_tube_49_box_membership_follows_stage_position(manifest):
    """The same regression one level down: the recorded box, not the condition."""
    assert boxes_by_section(manifest, 49) == {
        "01": "far_label",
        "04": "far_label",
        "02": "near_label",
        "03": "near_label",
    }


@pytest.mark.requires_data
def test_tube_45_regression_second_non_monotonic_slide(manifest):
    """TUBE 45 REGRESSION -- stage order ``03 04 | 01 02``.

    The second slide whose acquisition order does not follow slide position, and
    the more treacherous of the two: 45 is a ``both`` slide, so every section is
    positive and a wrong *condition* would not show up.  What can still go wrong
    is the recorded **box membership**, which feeds the level-matching logic
    (left-to-left anterior, right-to-right posterior) and any later
    re-derivation of condition.

    Near-label box (low stage X) = sections 03 and 04.  Far box = 01 and 02.
    """
    assert boxes_by_section(manifest, 45) == {
        "03": "near_label",
        "04": "near_label",
        "01": "far_label",
        "02": "far_label",
    }


@pytest.mark.requires_data
def test_tube_45_all_four_sections_positive(manifest):
    """45 is a ``both`` slide, so the non-monotonic layout must not flip anything."""
    assert conditions_by_section(manifest, 45) == {
        "01": "positive",
        "02": "positive",
        "03": "positive",
        "04": "positive",
    }


# ==========================================================================
# `both` slides -- animals with no negative control at all
# ==========================================================================
@pytest.mark.requires_data
@pytest.mark.parametrize("tube", sorted(BOTH_BOXES_TUBES))
def test_both_box_slides_have_no_negative_section(manifest, tube):
    """Tubes 35, 38, 45, 53 were stained in both boxes: 4+0 and 3+0.

    Marking any section of these animals ``negative`` would feed a fully-stained
    section into the negative-control gate, where it would fail the near-zero
    check and could trigger an exclusion -- of the positive data.
    """
    conditions = conditions_by_section(manifest, tube)
    assert set(conditions.values()) == {"positive"}, (
        f"tube {tube} has positive_box=both, so every section is positive; got "
        f"{conditions}"
    )


@pytest.mark.requires_data
@pytest.mark.parametrize("tube", sorted(BOTH_BOXES_TUBES))
def test_both_box_slides_report_has_negative_control_false(manifest, tube):
    """The flag the negative-control gate keys off must be False for these four.

    Spec section 2: the negative-control checks are defined *per
    animal-where-negatives-exist*.  Without this flag the gate either silently
    skips these animals or evaluates nothing and reports a pass.
    """
    values = {as_bool(v) for v in col(rows_for_tube(manifest, tube),
                                     "has_negative_control").tolist()}
    assert values == {False}, (
        f"tube {tube} has no negative control; has_negative_control={values}"
    )


@pytest.mark.requires_data
def test_animals_with_a_negative_box_report_has_negative_control_true(manifest):
    """The complement -- otherwise the flag could be hard-coded False."""
    for tube in (29, 30, 42, 49):
        values = {as_bool(v) for v in col(rows_for_tube(manifest, tube),
                                         "has_negative_control").tolist()}
        assert values == {True}, f"tube {tube} does have a negative box; got {values}"


@pytest.mark.requires_data
def test_exactly_four_animals_lack_a_negative_control(manifest):
    """A count check, so a future edit to slides.csv cannot widen this silently."""
    without = set()
    for tube in ALL_TUBES:
        rows = rows_for_tube(manifest, tube)
        values = {as_bool(v) for v in col(rows, "has_negative_control").tolist()}
        if values == {False}:
            without.add(tube)
    assert without == BOTH_BOXES_TUBES, (
        f"animals without a negative control changed: {sorted(without)} "
        f"(expected {sorted(BOTH_BOXES_TUBES)})"
    )


# ==========================================================================
# tube 37 -- a contradicted record must be excluded, never guessed
# ==========================================================================
@pytest.mark.requires_data
def test_tube_37_resolved_to_the_far_box(manifest):
    """TUBE 37 -- resolved at the bench on 2026-08-07.

    History worth keeping, because it is why the cross-check exists: slides.csv said
    far_label while the row's own annotation named section 01, which is in the NEAR
    box. The two sources named opposite boxes, and a coin flip there is a 50 % chance
    of measuring this animal's negative controls as data. The manifest refused to
    guess and marked every section ``unresolved`` until the bench answered.

    The bench answer: far_label is correct, positives are 03 and 04; the "use 01"
    note was a mis-write for "use 03". Independently confirmed from a QuPath
    thumbnail -- sections 01 and 02 render blue (DAPI only), 03 and 04 render
    colourful (marker signal present).
    """
    conditions = conditions_by_section(manifest, 37)
    assert len(conditions) == 4, f"tube 37 should contribute 4 rows: {conditions}"
    assert conditions["03"] == "positive"
    assert conditions["04"] == "positive"
    assert conditions["01"] == "negative"
    assert conditions["02"] == "negative"


@pytest.mark.requires_data
def test_unresolved_sections_are_reported_not_merely_excluded(manifest, manifest_mod):
    """The exclusion must be *reported*, not merely applied.

    A row quietly marked unresolved and never mentioned again is an animal that drops
    out of the analysis with no audit trail. Tube 37 was the live case and is now
    answered, so this asserts the behaviour against whatever is currently open in
    slides.csv -- and asserts silence when nothing is.
    """
    problems = manifest_mod.validate_manifest(manifest)
    assert isinstance(problems, list), (
        f"validate_manifest must return a list of messages, got {type(problems)}"
    )
    text = " ".join(str(p) for p in problems).lower()

    # Unresolved rows arise from two sources: a contradicted slides.csv row, and a
    # rescan whose sections cannot be matched to the original by stage X (tubes 33, 42
    # and 54, whose "rescan" holds a section the original does not have). Either way
    # the REQUIREMENT is the same and is what this asserts: refuse to guess, exclude
    # from the analysis manifest, and say so. Asserting a zero count instead would
    # turn a correct refusal into a red build.
    from _stage3_helpers import norm_condition
    unresolved_tubes = sorted({
        int(t) for t in col(manifest[col(manifest, "condition").map(norm_condition)
                                     == "unresolved"], "tube_id").unique()
    })
    if not unresolved_tubes:
        return

    for tube in unresolved_tubes:
        assert str(tube) in text, (
            f"tube {tube} has unresolved sections but validate_manifest said nothing "
            f"about it -- a silent drop-out. Messages: {problems}"
        )
    assert any(
        word in text for word in ("unresolved", "needs_confirmation", "confirm")
    ), f"the message does not say why it is a problem: {problems}"


@pytest.mark.requires_data
def test_validate_manifest_is_quiet_once_the_unresolved_rows_are_removed(
    manifest, manifest_mod
):
    """The complement: validate_manifest must not cry wolf on a clean manifest.

    A validator that always returns something trains its readers to ignore it,
    which is how the tube 37 message gets missed.
    """
    tube_column = resolve_column(manifest, "tube_id")
    clean = manifest[
        ~manifest[tube_column].astype("int64").isin(sorted(NEEDS_CONFIRMATION_TUBES))
    ]
    problems = manifest_mod.validate_manifest(clean)
    leftovers = [p for p in problems if "37" in str(p)]
    assert not leftovers, (
        f"tube 37 is no longer in the manifest but is still reported: {leftovers}"
    )


@pytest.mark.requires_data
def test_unresolved_rows_are_absent_now_the_bench_has_answered(manifest):
    """No other animal may be unresolved -- that would be a silent drop-out."""
    conditions = col(manifest, "condition")
    from _stage3_helpers import norm_condition

    unresolved = manifest[conditions.map(norm_condition) == "unresolved"]
    # The requirement is not "zero unresolved" -- new data legitimately raises new
    # questions. It is that every unresolved row is kept OUT of the measurement set.
    if len(unresolved):
        assert not unresolved["use_for_measurement"].any(), (
            "unresolved sections leaked into use_for_measurement: "
            f"{sorted(set(unresolved[unresolved['use_for_measurement']]['tube_id']))}"
        )


@pytest.mark.requires_data
def test_no_condition_value_is_unrecognised(manifest):
    """Every condition cell must be one of the three known states."""
    from _stage3_helpers import norm_condition

    bad = sorted({
        v for v in col(manifest, "condition").tolist()
        if norm_condition(v).startswith("?")
    }, key=repr)
    assert not bad, f"unrecognised condition values: {bad}"


# ==========================================================================
# condition must agree with slides.csv for every animal, not just the samples
# ==========================================================================
@pytest.mark.requires_data
def test_positive_section_count_matches_slides_csv_for_every_animal(manifest, slides):
    """Sweep: the number of positives per animal, derived independently.

    Computed here from ``slides.csv`` plus the box sizes the manifest itself
    recorded, so this catches a per-animal mistake anywhere in the cohort rather
    than only in the hand-picked cases above.
    """
    mismatches = {}
    for tube in ALL_TUBES:
        if tube in NEEDS_CONFIRMATION_TUBES:
            continue
        positive_box = slides[tube]["positive_box"]
        boxes = boxes_by_section(manifest, tube)
        conditions = conditions_by_section(manifest, tube)
        if positive_box == "both":
            expected = set(boxes)
        else:
            expected = {s for s, b in boxes.items() if b == positive_box}
        observed = {s for s, c in conditions.items() if c == "positive"}
        if expected != observed:
            mismatches[tube] = {
                "positive_box": positive_box,
                "boxes": boxes,
                "expected_positive": sorted(expected),
                "manifest_positive": sorted(observed),
            }
    assert not mismatches, f"condition disagrees with slides.csv: {mismatches}"


@pytest.mark.requires_data
def test_every_box_has_exactly_one_condition(manifest):
    """The rule that makes boxes meaningful: one box, one condition.

    Two sections in the same PAP-pen box cannot differ in staining -- they sat
    in the same well of antibody.  A manifest where they do has lost the box
    structure entirely.
    """
    from _stage3_helpers import norm_box, norm_condition

    offenders = {}
    for tube in ALL_TUBES:
        rows = rows_for_tube(manifest, tube)
        pairs = {}
        for box, condition in zip(
            col(rows, "box").map(norm_box), col(rows, "condition").map(norm_condition)
        ):
            pairs.setdefault(box, set()).add(condition)
        mixed = {b: c for b, c in pairs.items() if len(c) > 1}
        if mixed:
            offenders[tube] = mixed
    assert not offenders, f"a PAP-pen box holds two conditions: {offenders}"


# ==========================================================================
# index-only animals -- "metadata known, pixels not transferred" is normal
# ==========================================================================
@pytest.mark.requires_data
def test_index_only_animals_still_produce_rows(manifest):
    """23 of 31 animals have no payload folder. That is a state, not an error.

    If absent pixels raised or dropped the animal, the manifest would describe
    only the 8 animals that happen to be on this laptop, and the cohort would
    quietly shrink to a quarter of itself.
    """
    payload_column = find_column(manifest, "payload_present")
    assert payload_column is not None, (
        "the manifest has no payload_present column, so it cannot distinguish "
        "'pixels not transferred yet' from 'pixels missing'"
    )
    without = set()
    for tube in ALL_TUBES:
        rows = rows_for_tube(manifest, tube)
        assert len(rows) >= 3, f"tube {tube} contributed no rows"
        values = {as_bool(v) for v in rows[payload_column].tolist()}
        if values == {False}:
            without.add(tube)
    # Every tube contributing rows is the invariant, and it is asserted above.
    # Whether any animal is *currently* index-only is a fact about the transfer, not
    # about the code: it was 23 of 31, and as of 2026-08-08 it is none. Requiring at
    # least one would turn completing the cohort into a test failure.
    if not without:
        pytest.skip("all payloads have now been transferred; nothing index-only left")


@pytest.mark.requires_data
def test_payload_present_is_true_only_where_the_folder_exists(manifest, data_root):
    """payload_present must reflect the filesystem, not a hard-coded list."""
    payload_column = resolve_column(manifest, "payload_present")
    for tube in ALL_TUBES:
        rows = rows_for_tube(manifest, tube)
        values = {as_bool(v) for v in rows[payload_column].tolist()}
        on_disk = (Path(data_root) / f"_Image_{tube}_").is_dir()
        if not on_disk:
            assert values == {False}, (
                f"tube {tube} has no _Image_{tube}_ folder but the manifest says "
                f"payload_present={values}"
            )
    present = {t for t in ALL_TUBES if (Path(data_root) / f"_Image_{t}_").is_dir()}
    assert present <= set(PAYLOAD_TUBES), (
        f"unexpected payload folders on disk: {sorted(present - set(PAYLOAD_TUBES))}"
    )


@pytest.mark.requires_data
def test_a_missing_payload_does_not_blank_the_metadata(manifest):
    """Index-only rows still carry the metadata the index alone provides.

    The whole point of reading the ``.vsi`` index rather than going through
    Bio-Formats is that stage coordinates, pixel size and exposure survive
    without the payload.  A row that is present but empty is no better than a
    missing row.
    """
    payload_column = resolve_column(manifest, "payload_present")
    index_only = manifest[manifest[payload_column].map(as_bool) == False]  # noqa: E712
    if len(index_only) == 0:
        pytest.skip("all payloads have now been transferred; nothing index-only left")
    stage_column = find_column(manifest, "stage_x_um")
    if stage_column is not None:
        assert index_only[stage_column].notna().all(), (
            "index-only rows are missing stage X, so their box assignment cannot "
            "have been derived"
        )
    for canonical in ("section_label", "box", "condition"):
        column = resolve_column(manifest, canonical)
        assert index_only[column].notna().all(), (
            f"index-only rows have a null {canonical}"
        )


# ==========================================================================
# rescans
# ==========================================================================
@pytest.mark.requires_data
def test_rescans_are_flagged_and_excludable(manifest_mod, data_root, manifest):
    """Tubes 51 and 60 exist twice. Neither scan may impersonate the other.

    An unflagged rescan turns one animal into six rows for a four-section slide,
    which area-weighted aggregation would treat as extra tissue from the same
    mouse -- pseudo-replication that no downstream check would notice.
    """
    rescan_dir = Path(data_root) / "Rescan"
    if not rescan_dir.is_dir():
        pytest.skip("no Rescan/ directory under the raw-data root")

    without = manifest_mod.build_manifest(raw_root=data_root, include_rescans=False)
    assert len(without) < len(manifest), (
        "include_rescans=False produced the same number of rows as True, so the "
        "rescans of tubes 51 and 60 are either always included or never included"
    )

    mask = rescan_mask(manifest)
    assert mask is not None, (
        "rescan rows are present but nothing marks them. Either a boolean "
        "is_rescan column or a categorical scan column reading original/rescan "
        f"is needed. Columns: {list(manifest.columns)}"
    )
    flagged = manifest[mask]
    # Which tubes have a second scan is a fact about the DISK, not a constant. It began
    # as {51, 60} (re-acquired to fix exposure) and grew to {33, 42, 49, 51, 54, 60}
    # when spare booking time was used to image sections the first pass had missed.
    # Assert the property -- every flagged tube really does have a Rescan/ file -- not
    # a frozen membership list.
    import re as _re
    on_disk = {int(m.group(1)) for f in rescan_dir.glob("Image_*.vsi")
               if (m := _re.search(r"Image_(\d+)", f.name))}
    flagged_tubes = {int(t) for t in col(flagged, "tube_id").unique()}
    assert flagged_tubes <= on_disk, (
        f"manifest flags tubes {sorted(flagged_tubes - on_disk)} as rescanned, but "
        f"there is no file for them in {rescan_dir}"
    )
    assert len(flagged) > 0, "include_rescans=True added no flagged rows"
    assert len(manifest[~mask]) == len(without), (
        "dropping the flagged rows from the include_rescans=True manifest does "
        "not give the include_rescans=False manifest, so the flag does not "
        "actually identify the rescan rows"
    )


@pytest.mark.requires_data
def test_rescan_rows_do_not_disturb_the_original_rows(manifest_mod, data_root, manifest):
    """The original scan's rows must be identical with and without rescans."""
    rescan_dir = Path(data_root) / "Rescan"
    if not rescan_dir.is_dir():
        pytest.skip("no Rescan/ directory under the raw-data root")
    without = manifest_mod.build_manifest(raw_root=data_root, include_rescans=False)
    for tube in (51, 60):
        assert conditions_by_section(manifest, tube) == conditions_by_section(
            without, tube
        ), f"tube {tube}'s conditions changed when the rescan was folded in"


# ==========================================================================
# read-only raw data
# ==========================================================================
@pytest.mark.requires_data
def test_building_the_manifest_does_not_write_to_raw_data(manifest_mod, data_root):
    """Spec section 2: raw data is read-only. Ever.

    Cheap top-level check -- a new file, a deleted file or a touched index would
    all show up here.  It is not exhaustive, but the failure mode it guards
    against (a cache or a sidecar written next to the originals) is exactly the
    kind that lands at the top level.
    """
    root = Path(data_root)
    before = {p.name: p.stat().st_mtime for p in root.iterdir()}
    manifest_mod.build_manifest(raw_root=data_root, include_rescans=True)
    after = {p.name: p.stat().st_mtime for p in root.iterdir()}
    assert set(before) == set(after), (
        f"raw tree gained/lost entries: {set(after) ^ set(before)}"
    )
    touched = [n for n in before if before[n] != after[n]]
    assert not touched, f"build_manifest modified raw files: {touched}"


# ==========================================================================
# write_manifest
# ==========================================================================
@pytest.mark.requires_data
def test_write_manifest_writes_inside_out_dir_and_round_trips(
    manifest_mod, manifest, tmp_path
):
    """Artefacts land where they were asked to, and survive the round trip."""
    out_dir = tmp_path / "manifest_out"
    result = manifest_mod.write_manifest(manifest, out_dir)
    assert isinstance(result, dict), f"write_manifest returned {type(result)}"
    assert result, "write_manifest reported no artefacts"

    written = [Path(v) for v in result.values() if isinstance(v, (str, Path))]
    written = [p for p in written if p.suffix]
    assert written, f"write_manifest returned no file paths: {result}"
    for path in written:
        assert path.exists(), f"{path} was reported but does not exist"
        assert out_dir.resolve() in path.resolve().parents, (
            f"{path} was written outside the requested out_dir {out_dir}"
        )

    csvs = [p for p in written if p.suffix.lower() == ".csv"]
    assert csvs, f"no CSV among the written artefacts: {written}"
    reread = pd.read_csv(csvs[0])
    assert len(reread) == len(manifest), (
        f"wrote {len(manifest)} rows, read back {len(reread)}"
    )


# ==========================================================================
# crosscheck_condition_against_pixels -- auditor only, never a relabeller
# ==========================================================================
@pytest.mark.requires_data
@pytest.mark.requires_payload
@pytest.mark.slow
def test_crosscheck_never_relabels_a_section(manifest_mod, manifest, payload_datasets):
    """ADR-0003: pixels AUDIT slides.csv. They never overrule it.

    Positives and negatives are separable from the pixels with a large margin
    (GFAP 99.9th percentile 3358-5736 vs 196-388), which makes it tempting to
    "fix" a disagreement automatically.  Doing so would silently convert a
    failed positive stain into a fake negative control, and would make the
    negative-control check circular -- the controls would have been selected for
    the very property being tested.

    A disagreement is a blocking discrepancy for a human.  So: the input frame
    must come back untouched, and any condition the report carries must equal
    the one it was given.
    """
    before = manifest.copy(deep=True)
    report = manifest_mod.crosscheck_condition_against_pixels(
        manifest, level=3, max_slides=1
    )
    assert isinstance(report, pd.DataFrame), f"crosscheck returned {type(report)}"

    pd.testing.assert_frame_equal(
        manifest, before, check_like=True,
        obj="the manifest passed to crosscheck_condition_against_pixels",
    )

    condition_column = find_column(report, "condition")
    if condition_column is not None and len(report):
        tube_column = find_column(report, "tube_id")
        section_column = find_column(report, "section_label")
        if tube_column and section_column:
            from _stage3_helpers import norm_condition

            for _, row in report.iterrows():
                expected = conditions_by_section(manifest, int(row[tube_column]))
                key = section_key(row[section_column])
                if key in expected:
                    assert norm_condition(row[condition_column]) == expected[key], (
                        f"crosscheck changed tube {row[tube_column]} section {key} "
                        f"from {expected[key]} to {row[condition_column]!r}"
                    )

    assert set(report.columns) - set(manifest.columns), (
        "the crosscheck report adds no columns of its own, so it cannot be "
        "reporting a verdict"
    )


@pytest.mark.requires_data
@pytest.mark.requires_payload
@pytest.mark.slow
def test_crosscheck_respects_max_slides(manifest_mod, manifest, payload_datasets):
    """`max_slides` must bound the work -- it is the only way to run this cheaply."""
    report = manifest_mod.crosscheck_condition_against_pixels(
        manifest, level=3, max_slides=1
    )
    tube_column = find_column(report, "tube_id")
    if tube_column is None or not len(report):
        pytest.skip("crosscheck report carries no tube column to count")
    assert report[tube_column].nunique() <= 1, (
        f"max_slides=1 but the report covers {report[tube_column].nunique()} slides"
    )


# ==========================================================================
# environment guard
# ==========================================================================
def test_raw_data_root_is_not_written_by_these_tests():
    """These tests must never need write access to Dropbox."""
    assert os.environ.get("IHC_ALLOW_RAW_WRITES") is None, (
        "IHC_ALLOW_RAW_WRITES is set; raw data is read-only (spec section 2)"
    )
