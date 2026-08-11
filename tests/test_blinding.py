"""Tests for `ihc.ingest.blinding` -- these try to BREAK the blinding.

WHAT CLASS OF BUG THIS FILE PROTECTS AGAINST
============================================
**An artefact that is called blinded but is not.**

The delineator must stay blinded to treatment group. They
also did the staining and the imaging, so he can infer group from a tube ID
alone -- which means the blinded artefacts have to leak *nothing*, not merely
"nothing obvious".  Every leak vector below was verified present in this data:

* the file path and filename carry the tube ID (``Image_49.vsi``);
* ``stack1`` is the slide LABEL image -- it shows ``"1007344 - 29"`` as printed
  text **and** as a DataMatrix barcode encoding the same;
* the internal series names carry it too: tube 60's read
  ``"60_20x_DAPI, FITC, Cy3_01"``, and QuPath displays internal names in
  preference to filenames;
* the tube ID is in VSI tags 2061 and 120635;
* the acquisition timestamp is a *perfect* proxy -- scanning ran in ascending
  tube order, which is ascending group order;
* exposure identifies tubes 51 and 60, the only two acquired at non-standard
  settings.

And the structural trap: tube IDs run in contiguous treatment blocks --
29-40 Rapamycin Diet, 41-48 Extra Control Diet, 49-54 Control IP,
55-58 + 60 Rapamycin IP.  So **any order-preserving coding scheme reproduces
the group structure exactly**, including a sequential counter, a sort of the
file listing, or an arithmetic transform of the tube ID.  Codes must come from a
random permutation with a recorded seed.

A failed blinding does not announce itself.  Region boundaries drawn while
unblinded are still plausible boundaries; the numbers that follow are still in
range.  The damage only appears at peer review, if at all.

META-TESTS
==========
The tests below rest on three detectors: a substring/value scanner, a partition
comparison, and a block-structure statistic.  A broken detector would make every
blinding test pass *vacuously*, which is the worst possible outcome here.  So
each detector is itself tested against a deliberately-leaking synthetic frame
(``test_detector_*``), and those meta-tests run even when
``ihc.ingest.blinding`` does not exist yet.
"""

from __future__ import annotations

import hashlib
import os
import stat
from pathlib import Path

import pytest

from _stage3_helpers import (
    ALL_TUBES,
    GROUP_BLOCKS,
    WORK_DIR_FORBIDDEN,
    adjacency_null,
    associated_tube_id_hits,
    code_column,
    column_name_hits,
    exact_tube_id_hits,
    find_column,
    forbidden_substring_hits,
    group_of_tube,
    mean,
    partition_from,
    partition_of_groups,
    quantile,
    same_group_adjacent_pairs,
    stdev,
    substring_tube_id_hits,
    synthetic_private_manifest,
    text_forbidden_hits,
    text_tube_id_hits,
    tubes_of_blinded_rows,
)

pd = pytest.importorskip("pandas")

#: Seeds used wherever a test needs "several independent blindings".  Fixed, so
#: every probabilistic assertion below is deterministic: for a given
#: implementation the suite either always passes or always fails.
# Strong seeds, because production refuses weak ones (see _reject_guessable_seed) and
# these run through the same path. Any 12 distinct values do the job -- what these tests
# check is that different seeds give different mappings and that the block structure is
# destroyed, neither of which depends on the seed being memorable.
SEEDS = (
    100000000001, 100000000007, 100000000042, 100000000101,
    100000002718, 100000031415, 100000000813, 100000099991,
    100000123456, 100008675309, 100000000555, 100000000013,
)

#: The seed used for single-blinding tests.
SEED = SEEDS[0]


# ==========================================================================
# fixtures
# ==========================================================================
@pytest.fixture(scope="module")
def blinding():
    """`ihc.ingest.blinding`, or a clean skip while it is still being written."""
    return pytest.importorskip(
        "ihc.ingest.blinding",
        reason="ihc.ingest.blinding is not written yet (stage 3, in progress)",
    )


@pytest.fixture(scope="module")
def private_source():
    """A synthetic manifest carrying every known leak vector.

    Synthetic rather than real because blinding is pure table surgery -- reading
    31 index files would make the tests slower without making them stronger --
    and because a synthetic frame lets us plant the leaks deliberately.  A
    blinding test run against a manifest with nothing to leak proves nothing.
    """
    return synthetic_private_manifest(pd)


@pytest.fixture(scope="module")
def codes(blinding):
    return blinding.generate_codes(list(ALL_TUBES), seed=SEED)


@pytest.fixture(scope="module")
def split(blinding, private_source, codes):
    private_df, blinded_df = blinding.split_manifest(private_source, codes)
    assert isinstance(private_df, pd.DataFrame)
    assert isinstance(blinded_df, pd.DataFrame)
    return private_df, blinded_df


@pytest.fixture(scope="module")
def blinded(split):
    return split[1]


@pytest.fixture(scope="module")
def private(split):
    return split[0]


# ==========================================================================
# 1. the code mapping itself
# ==========================================================================
def test_same_seed_reproduces_the_same_mapping(blinding):
    """Reproducibility is the whole point of recording a seed.

    Without it the key cannot be regenerated if the custodian's file is lost,
    and the lock commit (spec section 10) cannot be verified after unblinding.
    """
    first = blinding.generate_codes(list(ALL_TUBES), seed=SEED)
    second = blinding.generate_codes(list(ALL_TUBES), seed=SEED)
    assert first == second


def test_the_same_seed_is_stable_under_input_order(blinding):
    """The mapping must depend on the seed, not on the order tubes arrive in.

    If it depends on iteration order, then a directory listing, a re-sorted
    slides.csv or a different filesystem silently produces a different key --
    and the recorded seed no longer reproduces the blinding.
    """
    ascending = blinding.generate_codes(list(ALL_TUBES), seed=SEED)
    shuffled = blinding.generate_codes(list(reversed(ALL_TUBES)), seed=SEED)
    assert ascending == shuffled, (
        "generate_codes gave a different mapping for the same tubes in a "
        "different order, so the seed alone does not determine the key"
    )


@pytest.mark.parametrize("other", SEEDS[1:5])
def test_a_different_seed_gives_a_different_mapping(blinding, other):
    """Otherwise the seed is decorative and the permutation is fixed."""
    if other == SEED:
        pytest.skip("same seed")
    assert blinding.generate_codes(list(ALL_TUBES), seed=other) != blinding.generate_codes(
        list(ALL_TUBES), seed=SEED
    )


def test_codes_are_a_bijection_over_the_cohort(blinding, codes):
    """One code per animal, all distinct, all strings."""
    assert set(codes) == set(ALL_TUBES), (
        f"codes cover {sorted(set(codes) ^ set(ALL_TUBES))} incorrectly"
    )
    values = list(codes.values())
    assert len(set(values)) == len(values), "two animals share a code"
    assert all(isinstance(v, str) for v in values), (
        f"codes must be strings; got {[type(v) for v in values[:3]]}"
    )


def test_a_code_never_contains_its_own_tube_id(codes):
    """`"A29"` is not a code, it is the tube ID with a letter in front."""
    offenders = {t: c for t, c in codes.items() if str(t) in str(c)}
    assert not offenders, f"codes containing their own tube ID: {offenders}"


def test_codes_are_not_a_hash_of_the_tube_id(codes):
    """A hash is deterministic without the key, so 31 guesses invert it.

    Spec section 2 bans "a hash or arithmetic transform of the tube ID"
    explicitly.  This checks the implementations somebody would actually reach
    for; the block-structure test below covers order-preserving schemes, which
    is the other half of the ban.
    """
    hits = []
    for tube, code in codes.items():
        text = str(code).lower()
        for algorithm in ("md5", "sha1", "sha256"):
            for payload in (str(tube), f"{tube:02d}", f"Image_{tube}"):
                digest = hashlib.new(algorithm, payload.encode()).hexdigest()
                if text and (digest.startswith(text) or text.startswith(digest[:6])):
                    hits.append((tube, code, algorithm, payload))
    assert not hits, f"codes look like hashes of the tube ID: {hits}"


# ==========================================================================
# 2. no tube ID anywhere in the blinded manifest
# ==========================================================================
def test_no_blinded_cell_equals_a_tube_id(blinded):
    """Exact-value scan across every cell of every column, all 31 IDs.

    Safe over numeric columns: no legitimate measurement in this dataset --
    stage coordinates in um, pixel size 0.325, widths around 16000, exposures in
    ms -- is exactly equal to an integer in 29..60.
    """
    hits = exact_tube_id_hits(blinded)
    assert not hits, (
        "the blinded manifest contains tube IDs as literal values:\n  "
        + "\n  ".join(f"column {c!r} row {r}: {v!r} == tube {t}" for c, r, v, t in hits[:20])
    )


def test_no_blinded_string_cell_contains_a_tube_id(blinded):
    """Substring scan of genuine string cells, anchored to digit boundaries.

    This is the one that catches ``"Image_49.vsi"``, ``"49_20x_DAPI..."`` and an
    ISO acquisition timestamp ending ``-07-29``.  A timestamp being flagged is a
    true positive: the cohort was scanned in ascending tube order, so the
    timestamp is a perfect group proxy and has no business in a blinded file.
    """
    hits = substring_tube_id_hits(blinded)
    assert not hits, (
        "the blinded manifest contains tube IDs inside string values:\n  "
        + "\n  ".join(f"column {c!r} row {r}: {v!r} contains {t}" for c, r, v, t in hits[:20])
    )


def test_no_blinded_cell_contains_a_group_name_or_a_path(blinded):
    """Group names, arm names, path separators, ``.vsi``, ``stack1``, ``60_``.

    ``60_`` is tube 60's internal series-name prefix.  QuPath shows internal
    series names, not filenames, so a project built from unrenamed series names
    displays the tube ID to the person who is meant to be blinded.
    """
    hits = forbidden_substring_hits(blinded)
    assert not hits, (
        "the blinded manifest contains unblinding tokens:\n  "
        + "\n  ".join(f"column {c!r} row {r}: {v!r} contains {t!r}" for c, r, v, t in hits[:20])
    )


def test_no_blinded_column_name_is_itself_a_leak(blinded):
    """A column called ``rapamycin_flag`` leaks before a single value is read."""
    hits = column_name_hits(blinded)
    assert not hits, f"blinded column names carry unblinding tokens: {hits}"


def test_the_named_leak_columns_are_gone(blinded):
    """The explicit strip list from `config/config.yaml`.

    ``blinding.strip_from_blinded_manifest``:
    tube_id, group, arm, original_path, acquisition_order.
    """
    for canonical in ("tube_id", "group", "arm"):
        assert find_column(blinded, canonical) is None, (
            f"the blinded manifest still has a {canonical!r} column"
        )
    leftovers = [
        c for c in blinded.columns
        if str(c).strip().lower() in {"original_path", "path", "file", "filename",
                                      "acquisition_order", "scan_order", "order"}
    ]
    assert not leftovers, f"blinded manifest still carries {leftovers}"


# ==========================================================================
# 3. no single column separates the groups
# ==========================================================================
def test_no_blinded_column_partition_equals_the_group_partition(blinded, codes):
    """Build each column's induced partition; none may equal the true groups.

    This is the test that catches a leak without having to name it -- a column
    that has been quietly relabelled (``batch``, ``cohort``, ``stain_run``,
    ``exposure_group``) is caught by its *structure* rather than its name.
    """
    truth = partition_of_groups([group_of_tube(t) for t in tubes_of_blinded_rows(blinded, codes)])
    offenders = []
    for column in blinded.columns:
        if partition_from(blinded[column].tolist()) == truth:
            offenders.append(str(column))
    assert not offenders, (
        f"these blinded columns partition the rows exactly like the treatment "
        f"groups do: {offenders}"
    )


def test_no_near_injective_blinded_column_reproduces_the_group_blocks(blinded, codes):
    """The stronger version, and the one that actually catches a timestamp.

    An exact-partition test is necessary but not sufficient: an acquisition
    timestamp is *unique per row*, so its partition is all singletons and it
    slips past.  What makes it a leak is that **sorting by it recovers the
    groups**, because scanning ran in ascending tube order.

    So: for every column with nearly one distinct value per animal -- the shape
    an order-preserving identifier has -- sort the animals by it and count
    adjacent same-group pairs.  Tube order gives 27 of 30.  A random permutation
    gives about 7.7.  Anything in the extreme upper tail of the Monte-Carlo null
    is an order-preserving proxy for treatment group.
    """
    tubes = tubes_of_blinded_rows(blinded, codes)
    frame = blinded.copy()
    frame["_tube"] = tubes
    per_animal = frame.groupby("_tube", sort=True).first().reset_index()
    animal_tubes = per_animal["_tube"].tolist()
    groups = [group_of_tube(t) for t in animal_tubes]
    null = adjacency_null(groups)
    ceiling = quantile(null, 0.999)

    offenders = {}
    for column in per_animal.columns:
        if column == "_tube":
            continue
        values = per_animal[column].tolist()
        if len(set(map(repr, values))) < 0.8 * len(values):
            continue  # too many ties to be an ordering; the partition test covers it
        order = sorted(range(len(values)), key=lambda i: (_sort_key(values[i]), i))
        observed = same_group_adjacent_pairs([groups[i] for i in order])
        if observed > ceiling:
            offenders[str(column)] = (observed, ceiling)
    assert not offenders, (
        "sorting the animals by these blinded columns reproduces the treatment "
        f"block structure (observed vs 99.9th-percentile null): {offenders}"
    )


def _sort_key(value):
    """Total order over mixed-type cells, numbers before strings."""
    try:
        return (0, float(value), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(value))


# ==========================================================================
# 4. the codes destroy the block structure
# ==========================================================================
def _block_statistic(codes: dict) -> int:
    """Same-group adjacent pairs when animals are listed in code order."""
    ordered = sorted(codes, key=lambda tube: str(codes[tube]))
    return same_group_adjacent_pairs([group_of_tube(t) for t in ordered])


def test_code_order_does_not_reproduce_the_treatment_blocks(blinding):
    """The headline structural test.

    Tube IDs run in contiguous blocks of 12, 8, 6 and 5.  Listing the animals in
    tube order gives 27 same-group adjacencies out of 30; a random permutation
    averages 7.7.  Any order-preserving code -- sequential, file-order,
    ``tube - 28``, ``tube * 7 mod 61`` -- lands near 27 and hands the block
    structure straight back.

    Checked across twelve seeds so that one unlucky-but-legitimate permutation
    cannot fail the suite, while a systematically order-preserving scheme fails
    on all twelve.
    """
    groups = [group_of_tube(t) for t in ALL_TUBES]
    null = adjacency_null(groups)
    ceiling = quantile(null, 0.999)
    null_mean, null_sd = mean(null), stdev(null)

    observed = [
        _block_statistic(blinding.generate_codes(list(ALL_TUBES), seed=s)) for s in SEEDS
    ]
    extreme = [(s, o) for s, o in zip(SEEDS, observed) if o > ceiling]
    assert len(extreme) <= 1, (
        f"{len(extreme)} of {len(SEEDS)} seeds put the code order in the extreme "
        f"upper tail of the null (>{ceiling}); the coding scheme preserves the "
        f"treatment block structure. Offending seeds: {extreme}"
    )

    tolerance = 4 * null_sd / (len(SEEDS) ** 0.5)
    assert abs(mean(observed) - null_mean) <= max(tolerance, 1.0), (
        f"mean same-group adjacency over {len(SEEDS)} seeds is {mean(observed):.2f}, "
        f"against a null mean of {null_mean:.2f} (sd {null_sd:.2f}). Codes are not "
        f"behaving like a random permutation."
    )


def test_code_order_is_not_tube_order(blinding):
    """The degenerate cases, stated plainly so the failure message is obvious."""
    for seed in SEEDS[:4]:
        mapping = blinding.generate_codes(list(ALL_TUBES), seed=seed)
        by_code = sorted(mapping, key=lambda t: str(mapping[t]))
        assert by_code != list(ALL_TUBES), f"seed {seed}: codes are in tube order"
        assert by_code != list(reversed(ALL_TUBES)), (
            f"seed {seed}: codes are in reverse tube order"
        )


# ==========================================================================
# 5. the private manifest is the complement -- it MUST unblind
# ==========================================================================
def test_the_private_manifest_keeps_tube_id_and_group(private):
    """Without this the key is useless and nobody can ever unblind.

    The complement test matters as much as the leak tests: an implementation
    that strips everything from both frames passes every test above.
    """
    tube_column = find_column(private, "tube_id")
    group_column = find_column(private, "group")
    assert tube_column is not None, "the private manifest has no tube_id column"
    assert group_column is not None, "the private manifest has no group column"
    assert set(int(t) for t in private[tube_column].unique()) == set(ALL_TUBES)
    assert set(private[group_column].unique()) == set(GROUP_BLOCKS)


def test_the_private_manifest_carries_the_code(private, codes):
    """It is the join between the blinded world and the real one."""
    column = code_column(private, codes)
    tube_column = find_column(private, "tube_id")
    for tube, code in zip(private[tube_column].tolist(), private[column].tolist()):
        assert str(code) == str(codes[int(tube)]), (
            f"private manifest maps tube {tube} to {code!r}, key says {codes[int(tube)]!r}"
        )


def test_both_frames_describe_the_same_rows(private, blinded, private_source):
    """Splitting must not drop or duplicate a section."""
    assert len(private) == len(private_source)
    assert len(blinded) == len(private_source)


# ==========================================================================
# 6. unresolved sections must not be blinded
# ==========================================================================
def test_blinding_a_manifest_with_unresolved_sections_raises(blinding, codes):
    """Tube 37's condition is contradicted, so it has no place downstream.

    If an unresolved section reaches a blinded QuPath project it gets
    delineated, classified and measured, and the only record that it should not
    have been is a column that was stripped on the way in.  Refuse loudly at the
    boundary instead.
    """
    tainted = synthetic_private_manifest(
        pd,
        condition_of=lambda tube, label, condition: (
            "unresolved" if tube == 37 else condition
        ),
    )
    with pytest.raises((ValueError, RuntimeError)) as excinfo:
        blinding.split_manifest(tainted, codes)
    message = str(excinfo.value).lower()
    assert "unresolved" in message or "37" in message, (
        f"the refusal does not say what is wrong: {excinfo.value!r}"
    )


def test_a_clean_manifest_still_blinds(blinding, codes, private_source):
    """The complement -- otherwise split_manifest could just always raise."""
    private_df, blinded_df = blinding.split_manifest(private_source, codes)
    assert len(blinded_df) == len(private_source)


# ==========================================================================
# 7. audit_blinded
# ==========================================================================
def _leaks(problems):
    """Only the findings that mean the blinding is broken.

    `audit_blinded` returns two severities, `LEAK:` and `RISK:`.  A RISK on clean
    output is expected and is the honest answer, not a bug: `exposure_is_standard`
    really does single out the two animals the spec names in writing, and no
    amount of column-stripping changes that -- it is a property of the cohort.
    Asserting `problems == []` would force that finding to be suppressed, which
    is the opposite of what anyone wants.
    """
    return [p for p in problems if str(p).strip().upper().startswith("LEAK")]


def test_audit_passes_a_clean_blinded_manifest(blinding, blinded, codes, private):
    problems = blinding.audit_blinded(blinded, codes, private)
    assert isinstance(problems, list), f"audit_blinded returned {type(problems)}"
    assert _leaks(problems) == [], (
        f"audit_blinded reported a LEAK in its own output: {_leaks(problems)}"
    )
    for problem in problems:
        assert str(problem).strip().upper().startswith(("LEAK", "RISK")), (
            f"audit finding without a severity prefix: {problem!r}"
        )


def test_audit_catches_a_planted_tube_id(blinding, blinded, codes, private):
    """The auditor must fail when the blinding fails -- otherwise it is scenery."""
    tampered = blinded.copy()
    tampered["scan_note"] = [f"slide {t}" for t in tubes_of_blinded_rows(blinded, codes)]
    problems = blinding.audit_blinded(tampered, codes, private)
    assert _leaks(problems), (
        f"a column full of tube IDs was not reported as a LEAK: {problems}"
    )


def test_audit_catches_a_planted_group_column(blinding, blinded, codes, private):
    tampered = blinded.copy()
    tampered["batch"] = [group_of_tube(t) for t in tubes_of_blinded_rows(blinded, codes)]
    problems = blinding.audit_blinded(tampered, codes, private)
    assert _leaks(problems), (
        f"a column holding the treatment group was not reported as a LEAK: {problems}"
    )


def test_audit_catches_an_acquisition_timestamp(blinding, blinded, codes, private):
    """The leak an exact-partition test cannot see.

    Timestamps are unique per row, so the partition they induce is all
    singletons -- never equal to the group partition.  What makes them fatal is
    that they SORT into tube order, because scanning ran in ascending tube order
    and tube IDs run in contiguous treatment blocks.
    """
    tampered = blinded.copy()
    tubes = tubes_of_blinded_rows(blinded, codes)
    tampered["scanned_at"] = [
        f"2026-07-{20 + t // 24:02d}T{t % 24:02d}:00:00+00:00" for t in tubes
    ]
    problems = blinding.audit_blinded(tampered, codes, private)
    assert _leaks(problems), (
        f"a monotone acquisition timestamp was not reported as a LEAK: {problems}"
    )


# ==========================================================================
# 8. write_blinded -- where the leak becomes a file on disk
# ==========================================================================
@pytest.fixture
def written(blinding, private, blinded, codes, tmp_path):
    """Run write_blinded into two fresh directories and return them."""
    custodian_dir = tmp_path / "custodian"
    work_dir = tmp_path / "work"
    result = blinding.write_blinded(
        private, blinded, codes, SEED, custodian_dir=custodian_dir, work_dir=work_dir
    )
    assert isinstance(result, dict), f"write_blinded returned {type(result)}"
    return custodian_dir, work_dir, result


def _text_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        try:
            yield path, path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue


def test_nothing_in_the_work_directory_leaks_a_tube_id(written, codes):
    """The work tree is shared storage. Every byte of it is visible to the delineator.

    Scanned for *associated* tube IDs rather than bare digits.  A legitimate
    report says ``"n_animals": 31`` and ``"29 of 31 slides"`` -- the cohort size
    and a fact from the published spec, neither of which names a mouse.  A
    blanket digit scan flags both, gets muted within a week, and then misses the
    real thing.  What counts as a leak is a tube ID attached to something: a
    code on the same line, the word tube/animal/slide, or an ``Image_NN``
    filename stem.
    """
    _, work_dir, _ = written
    offenders = {}
    for path, text in _text_files(work_dir):
        hits = associated_tube_id_hits(text, codes)
        if hits:
            offenders[str(path.relative_to(work_dir))] = hits
    assert not offenders, f"tube IDs found in the shared work tree: {offenders}"


def test_nothing_in_the_work_directory_leaks_a_group_or_a_path(written):
    _, work_dir, _ = written
    offenders = {}
    for path, text in _text_files(work_dir):
        hits = text_forbidden_hits(text, WORK_DIR_FORBIDDEN)
        hits += text_forbidden_hits(text, ("rapamycin", "vehicle", "extra control"))
        if hits:
            offenders[str(path.relative_to(work_dir))] = sorted(set(hits))
    assert not offenders, f"unblinding tokens in the shared work tree: {offenders}"


def test_the_seed_is_never_written_to_the_work_directory(written):
    """The seed regenerates the whole permutation, so it IS the key.

    Recording it next to the blinded manifest -- in a provenance header, a
    "parameters used" block, a log line -- undoes the blinding for anyone who
    can run ``generate_codes``.
    """
    _, work_dir, _ = written
    offenders = [
        str(path.relative_to(work_dir))
        for path, text in _text_files(work_dir)
        if str(SEED) in text
    ]
    assert not offenders, f"the blinding seed appears in the shared work tree: {offenders}"


def test_the_custodian_directory_holds_the_key(written, codes):
    """The complement: the key has to exist somewhere, or nobody can unblind."""
    custodian_dir, _, _ = written
    assert custodian_dir.is_dir(), "write_blinded did not create the custodian directory"
    found = False
    for _, text in _text_files(custodian_dir):
        if any(str(t) in text for t in ALL_TUBES) and any(
            str(c) in text for c in codes.values()
        ):
            found = True
            break
    assert found, (
        "no file under the custodian directory contains both tube IDs and codes, "
        "so the mapping has not been recorded anywhere"
    )


def test_custodian_and_work_trees_are_disjoint(written):
    """The key must not sit inside the tree the blinded analyst works from."""
    custodian_dir, work_dir, _ = written
    assert custodian_dir.resolve() not in work_dir.resolve().parents
    assert work_dir.resolve() not in custodian_dir.resolve().parents
    assert custodian_dir.resolve() != work_dir.resolve()


def test_the_custodian_directory_is_not_world_readable(written):
    """ADR-0001 / `config/paths.yaml`: the custodian tree is mode 700.

    It holds the key on local disk, outside git so it cannot be pushed and
    outside Dropbox so it cannot be shared by accident. Group and other
    permissions defeat the point.
    """
    custodian_dir, _, _ = written
    mode = stat.S_IMODE(os.stat(custodian_dir).st_mode)
    assert mode & 0o077 == 0, (
        f"custodian directory is mode {mode:o}; it must be 700 (no group, no other)"
    )


def test_write_blinded_reports_where_it_put_things(written):
    """The returned dict has to be usable -- every path in it must exist."""
    custodian_dir, work_dir, result = written
    paths = [Path(v) for v in result.values() if isinstance(v, (str, Path))]
    paths = [p for p in paths if p.suffix or p.is_dir()]
    assert paths, f"write_blinded returned no paths: {result}"
    for path in paths:
        assert path.exists(), f"{path} was reported but does not exist"
        roots = [custodian_dir.resolve(), work_dir.resolve()]
        assert any(r == path.resolve() or r in path.resolve().parents for r in roots), (
            f"{path} is outside both the custodian and work directories"
        )


# ==========================================================================
# 9. integration against the real manifest, when it exists
# ==========================================================================
@pytest.mark.requires_data
def test_the_real_manifest_blinds_without_leaking(blinding, data_root):
    """End to end on the actual 31 index files, once stage 3 is complete.

    The synthetic tests above are sharper, but they only contain the leaks
    somebody thought to plant.  This one runs the same scanners over whatever
    the real pipeline actually produces.
    """
    manifest_mod = pytest.importorskip(
        "ihc.ingest.manifest", reason="ihc.ingest.manifest is not written yet"
    )
    from _stage3_helpers import (
        NEEDS_CONFIRMATION_TUBES,
        require_cohort_index_files,
        resolve_column,
    )

    require_cohort_index_files(data_root)
    df = manifest_mod.build_manifest(raw_root=data_root, include_rescans=True)
    tube_column = resolve_column(df, "tube_id")
    usable = df[~df[tube_column].astype("int64").isin(sorted(NEEDS_CONFIRMATION_TUBES))]

    tubes = sorted(set(int(t) for t in usable[tube_column].unique()))
    mapping = blinding.generate_codes(tubes, seed=SEED)
    _, blinded_df = blinding.split_manifest(usable, mapping)

    assert not exact_tube_id_hits(blinded_df, tubes)
    assert not substring_tube_id_hits(blinded_df, tubes)
    assert not forbidden_substring_hits(blinded_df)
    assert not column_name_hits(blinded_df)

    truth = partition_of_groups(
        [group_of_tube(t) for t in tubes_of_blinded_rows(blinded_df, mapping)]
    )
    for column in blinded_df.columns:
        assert partition_from(blinded_df[column].tolist()) != truth, (
            f"column {column!r} of the real blinded manifest separates the groups"
        )


# ==========================================================================
# META-TESTS: the detectors must detect
# --------------------------------------------------------------------------
# These run whether or not `ihc.ingest.blinding` exists.  Every assertion above
# is of the form "the scanner found nothing"; if a scanner can never find
# anything, all of them pass while the blinding is wide open.
# ==========================================================================
def _leaky_frame():
    """A blinded-looking frame that leaks in four different ways."""
    tubes = list(ALL_TUBES)
    return pd.DataFrame(
        {
            "coded_id": [f"XK{i:03d}" for i in range(len(tubes))],
            "tube_id": tubes,                                        # exact value
            "source": [f"/RawData/Image_{t}.vsi" for t in tubes],    # substring + path
            "cohort": [group_of_tube(t) for t in tubes],             # group partition
            "scanned_at": [f"2026-07-29T10:{i:02d}:00" for i in range(len(tubes))],
        }
    )


def test_detector_exact_value_scan_catches_a_planted_tube_id_column():
    hits = exact_tube_id_hits(_leaky_frame())
    assert hits, "the exact-value scanner missed a whole column of tube IDs"
    assert {c for c, _, _, _ in hits} >= {"tube_id"}


def test_detector_substring_scan_catches_a_planted_path():
    hits = substring_tube_id_hits(_leaky_frame())
    assert {c for c, _, _, _ in hits} >= {"source"}, (
        "the substring scanner missed 'Image_49.vsi'"
    )


def test_detector_substring_scan_does_not_fire_on_legitimate_numbers():
    """The false-positive guard, without which the scanner gets switched off.

    ``0.325`` contains "32", ``16384`` contains "38", ``1840.0`` contains "40".
    None of them is a tube ID.
    """
    clean = pd.DataFrame(
        {
            "coded_id": ["QZ118", "SLK2291", "MM7734"],
            "pixel_size_um": [0.325, 0.32501, 0.325],
            "width_px": [16384, 15360, 22528],
            "exposure_cy3_ms": [1840.0, 1840.0, 397.93],
            "stage_x_um": [12934.5, 29001.25, 45600.0],
        }
    )
    assert not exact_tube_id_hits(clean), "exact scanner fired on legitimate numbers"
    assert not substring_tube_id_hits(clean), (
        "substring scanner fired on legitimate numeric values -- it would be "
        "disabled within a week and the real leaks would go unnoticed"
    )


def test_detector_forbidden_token_scan_catches_group_names():
    hits = forbidden_substring_hits(_leaky_frame())
    assert {c for c, _, _, _ in hits} >= {"cohort", "source"}


def test_detector_partition_comparison_catches_a_relabelled_group_column():
    frame = _leaky_frame()
    truth = partition_of_groups([group_of_tube(t) for t in ALL_TUBES])
    assert partition_from(frame["cohort"].tolist()) == truth, (
        "the partition comparison cannot even see a verbatim group column"
    )
    assert partition_from(frame["coded_id"].tolist()) != truth


def test_detector_block_statistic_separates_random_from_sequential():
    """The statistic must put a sequential coder far above a random one.

    Sequential coding of contiguous treatment blocks gives 27 of a possible 30
    same-group adjacencies.  If this assertion ever fails, the block-structure
    test above is measuring nothing.
    """
    groups = [group_of_tube(t) for t in ALL_TUBES]
    null = adjacency_null(groups)
    ceiling = quantile(null, 0.999)

    sequential = {t: f"C{i:03d}" for i, t in enumerate(ALL_TUBES)}
    assert _block_statistic(sequential) == 27
    assert _block_statistic(sequential) > ceiling, (
        f"a sequential coder scored {_block_statistic(sequential)} against a "
        f"99.9th-percentile null of {ceiling}; the statistic has no power"
    )

    arithmetic = {t: f"C{(t * 7) % 61:03d}" for t in ALL_TUBES}
    assert _block_statistic(arithmetic) >= 0  # sanity: it is computable

    import random as _random

    rng = _random.Random(4242)
    shuffled = list(ALL_TUBES)
    rng.shuffle(shuffled)
    random_codes = {t: f"C{i:03d}" for i, t in enumerate(shuffled)}
    assert _block_statistic(random_codes) <= ceiling, (
        "a genuinely random permutation was flagged; the test would be flaky"
    )


def test_detector_text_scanners_catch_a_leaky_file(tmp_path):
    path = tmp_path / "blinded.csv"
    path.write_text("coded_id,note\nXK001,from Image_49.vsi (Rapamycin Diet)\n")
    text = path.read_text()
    assert text_tube_id_hits(text) == ["49"]
    assert set(text_forbidden_hits(text, WORK_DIR_FORBIDDEN)) >= {".vsi", "image_"}
    assert "rapamycin" in text_forbidden_hits(text, ("rapamycin",))


def test_detector_association_scan_catches_a_key_line(tmp_path):
    """A line holding a code and a tube ID together IS the key."""
    codes = {49: "XK007", 29: "QZ118"}
    text = "code,tube\nXK007,49\n"
    assert associated_tube_id_hits(text, codes), (
        "a line mapping a code to a tube number was not detected"
    )
    assert associated_tube_id_hits("slide 49 was rescanned", codes), (
        "a tube ID labelled with the word 'slide' was not detected"
    )
    assert associated_tube_id_hits("see _Image_49_/stack10002", codes), (
        "an Image_NN filename stem was not detected"
    )


def test_detector_association_scan_ignores_counts_and_prose(tmp_path):
    """The false-positive guard that keeps the scan switched on.

    Every string here is in the real blinding report, and none of them names an
    animal: 31 is the cohort size, 29 of 31 is the exposure split from the
    published spec, and 118 is a row count.
    """
    codes = {49: "XK007", 29: "QZ118"}
    for benign in (
        '"n_animals": 31',
        '"n_rows": 118',
        "a 29/1/1 split is not anonymous",
        '"adjacent_same_group_null_mean": 7.696',
        '"pixel_size_um": 0.325',
        '"width_px": 16384',
    ):
        assert associated_tube_id_hits(benign, codes) == [], (
            f"the association scan fired on {benign!r}, which names no animal"
        )


def test_detector_text_scanners_are_quiet_on_a_clean_file(tmp_path):
    path = tmp_path / "clean.csv"
    path.write_text(
        "coded_id,section_label,condition,pixel_size_um,width_px\n"
        "XK001,01,positive,0.325,16384\n"
        "SLK2291,02,negative,0.32501,15360\n"
    )
    text = path.read_text()
    assert text_tube_id_hits(text) == []
    assert text_forbidden_hits(text, WORK_DIR_FORBIDDEN) == []


def test_the_cohort_constants_still_describe_the_cohort():
    """Guards the assumption every statistic above rests on."""
    assert len(ALL_TUBES) == 31
    assert 59 not in ALL_TUBES
    sizes = sorted(len(v) for v in GROUP_BLOCKS.values())
    assert sizes == [5, 6, 8, 12]
    assert sum(sizes) == 31
    for tubes in GROUP_BLOCKS.values():
        assert all(isinstance(t, int) for t in tubes)


# ==========================================================================
# The seed guard
# ==========================================================================
# This section exists because the guard was DOCUMENTED before it was written. An ADR
# and a commit message both stated that weak seeds were rejected while the code had no
# such check, and nobody noticed for a week — a false assurance is worse than a known
# gap, because a documented control does not get re-checked. These tests make the claim
# and the code fail together.

@pytest.mark.parametrize("seed", [20260807, 19991231, 20000101])
def test_date_shaped_seeds_are_refused(blinding, seed):
    """A date is ~365 guesses.

    The mapping is deterministic in (sorted cohort, seed); the roster is a public
    constant in this repository and the issued codes are in the blinded manifest an
    analyst holds. So an attacker enumerates seeds and compares — no cryptanalysis
    required.
    """
    with pytest.raises(ValueError, match="date-shaped"):
        blinding.generate_codes(list(ALL_TUBES), seed=seed)


@pytest.mark.parametrize("seed", [1, 42, 999, 12345, 99999999999])
def test_short_numeric_seeds_are_refused(blinding, seed):
    """Anything enumerable in seconds is not a key."""
    with pytest.raises(ValueError, match="digits"):
        blinding.generate_codes(list(ALL_TUBES), seed=seed)


@pytest.mark.parametrize("seed", ["password", "secret", "blinding", "abc123"])
def test_weak_word_seeds_are_refused(blinding, seed):
    """A wordlist entry is not a key either."""
    with pytest.raises(ValueError, match="too simple"):
        blinding.generate_codes(list(ALL_TUBES), seed=seed)


def test_a_proper_random_seed_is_accepted(blinding):
    """secrets.randbits(128) is the documented way to draw one, and must work."""
    import secrets
    codes = blinding.generate_codes(list(ALL_TUBES), seed=secrets.randbits(128))
    assert len(codes) == len(ALL_TUBES)


def test_the_escape_hatch_is_opt_in_and_never_used_in_production(blinding):
    """Tests may bypass the guard; the shipped code may not.

    A weak seed is fine in a test that needs a readable constant. It is not fine
    anywhere that writes a real key, so the CLI must never pass the flag.
    """
    assert blinding.generate_codes(list(ALL_TUBES), seed=42, allow_weak_seed=True)
    cli = (Path(__file__).resolve().parents[1] / "ihc").read_text()
    assert "allow_weak_seed" not in cli, (
        "the ihc entry command passes allow_weak_seed — production must not bypass "
        "the seed guard"
    )


def test_timestamps_do_not_read_as_a_tube_id_leak():
    """A timestamp is not a leak, however many tube numbers it happens to contain.

    `"created_utc": "2026-08-11T15:34:49Z"` carries the tube tokens 34 and 49, and
    the substring T15 collides with the coded-ID label pool -- enough for the scanner
    to report "a tube ID on the same line as a code". It did, and it did so as a
    function of the time of day the artefact was written. An intermittent leak test
    is worse than none, because it gets muted.
    """
    codes = {29: "T15", 30: "B16"}
    assert associated_tube_id_hits('"created_utc": "2026-08-11T15:34:49Z",', codes) == []
    assert associated_tube_id_hits('"generated": "2026-08-11 15:34:49",', codes) == []
    # ... but a real association on the same line must still be caught.
    assert associated_tube_id_hits('"T15": 49,', codes)
    assert associated_tube_id_hits('tube 49 -> T15', codes)
