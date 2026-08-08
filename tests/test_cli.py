"""Smoke tests for the ``./ihc`` entry command.

Protects against: a subcommand that is completely dead while the unit tests all pass.

This file exists because exactly that happened. ``./ihc meta`` shipped calling
``assign_boxes(meta)`` where the function takes ``assign_boxes(meta.series)``, so it
raised ``TypeError`` on all 31 slides and reported "0 slides swept, 31 failed" — while
453 unit tests passed, because nothing exercised the CLI. These tests are deliberately
shallow: they check each subcommand runs, exits sensibly, and emits the one line that
proves it did its job. Depth belongs in the module tests.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
IHC = REPO / "ihc"


def run(*args, timeout=300):
    """Run ./ihc with *args and return the CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(IHC), *args],
        capture_output=True, text=True, timeout=timeout, cwd=REPO,
    )


def test_entry_command_exists_and_is_executable():
    assert IHC.exists(), "./ihc entry command is missing"


def test_no_args_prints_usage_and_does_not_crash():
    r = run()
    assert r.returncode in (0, 1, 2), f"unexpected exit {r.returncode}\n{r.stderr}"
    assert "ihc" in (r.stdout + r.stderr).lower()


@pytest.mark.parametrize("sub", ["doctor", "check-config", "meta", "verify"])
def test_subcommand_is_registered(sub):
    """A registered subcommand must not die with a Python traceback."""
    r = run(sub, "--help")
    combined = r.stdout + r.stderr
    assert "Traceback" not in combined, f"./ihc {sub} --help raised:\n{combined}"


def test_doctor_reports_interpreter_and_paths():
    r = run("doctor")
    assert "Traceback" not in r.stdout + r.stderr
    assert "interpreter" in r.stdout
    assert "paths" in r.stdout


def test_check_config_parses_all_three_yaml_files():
    r = run("check-config")
    out = r.stdout
    assert "Traceback" not in out + r.stderr
    for name in ("paths.yaml", "config.yaml", "channels.yaml"):
        assert f"{name} parses" in out, f"check-config did not confirm {name}"


def test_check_config_surfaces_rows_needing_confirmation():
    """Any row with an open question must be shouted about, never allowed to go quiet.

    Tubes 49 and then 37 were the live cases; both are now answered, so the set is
    empty and there is nothing to surface. The test reads slides.csv and asserts the
    surfacing behaviour only when there is something to surface — a test that
    demanded a warning forever would fail the moment the bench did its job.
    """
    import csv
    from pathlib import Path
    rows = list(csv.DictReader(open(Path(REPO) / "config" / "slides.csv")))
    open_rows = [r for r in rows if (r.get("needs_confirmation") or "").strip()]
    r = run("check-config")
    if open_rows:
        assert "NEEDS CONFIRMATION" in r.stdout or "needs_confirmation" in r.stdout, \
            f"{len(open_rows)} row(s) need confirmation but check-config was silent"
    else:
        assert "Traceback" not in r.stdout + r.stderr


@pytest.mark.slow
def test_meta_sweeps_every_slide_and_writes_an_artefact(data_root):
    """The regression test for the dead-subcommand bug.

    ``meta`` derives PAP-pen box membership, which is the science-critical rule in the
    whole spec. It must sweep every index file, fail none, and leave a file behind —
    transient terminal output is not an artefact anyone can diff or review.
    """
    if data_root is None or not Path(data_root).exists():
        pytest.skip("raw data root not available")
    r = run("meta", timeout=600)
    out = r.stdout
    assert "Traceback" not in out + r.stderr, out + r.stderr
    assert "0 failed" in out, f"meta reported failures:\n{out[-2000:]}"
    assert "slides swept" in out
    assert "series_metadata.csv" in out, "meta produced no artefact"
    # the tube 49 regression, end to end through the CLI
    assert "02,03" in out and "01,04" in out, \
        "tube 49 box assignment (02,03 | 01,04) not present in the meta output"
