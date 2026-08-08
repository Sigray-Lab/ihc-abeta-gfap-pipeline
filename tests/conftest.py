"""Shared fixtures and path resolution for the IHC pipeline test suite.

This file exists to keep two classes of bug out of the tests themselves:

1.  **Silent data-root drift.**  Tests that hard-code a Dropbox path pass on the
    maintainer's laptop and are meaningless anywhere else.  `data_root` resolves
    the raw-data location from an explicit precedence chain and *skips loudly*
    when it cannot, so an absent data root can never masquerade as a pass.

2.  **Dehydrated Dropbox placeholders.**  A cloud-only file reports a plausible
    `st_size` from `stat()` but stalls or fails when read.  Any fixture that
    hands out a payload folder checks `st_blocks` against `st_size` first and
    skips with an explicit message, rather than letting a test hang or fail for
    a reason that has nothing to do with the code under test.

Data-root precedence (first hit wins):
    1. ``$IHC_DATA_ROOT``
    2. ``config/paths.yaml`` in the repo root (keys tried in order:
       ``raw_data``, ``raw_data_root``, ``data_root``, ``raw``, ``rawdata``,
       optionally nested under a top-level ``paths:`` mapping)
    3. ``DEFAULT_DATA_ROOT`` below -- the documented project location
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# --------------------------------------------------------------------------
# make `src/ihc/...` importable without requiring an editable install
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


# --------------------------------------------------------------------------
# data-root resolution
# --------------------------------------------------------------------------
DEFAULT_DATA_ROOT = Path(
    "~/Dropbox/longevity_project/projects/Rapamycin_CSF_mice_Per"
    "/IHC_analysis_pipeline/RawData"
)

_PATHS_YAML = REPO_ROOT / "config" / "paths.yaml"

# Where the raw-data root lives in config/paths.yaml.  `roots.raw_root` is the
# current layout; the rest are accepted so a rename does not silently drop the
# whole suite back to the hard-coded fallback.
_YAML_ROOT_KEYS = ("raw_root", "raw_data", "raw_data_root", "data_root", "rawdata")

# The full expected index-file set: tubes 29-58 plus 60.  Mouse 59 was excluded
# before imaging (mounting fault) and has no .vsi at all.
ALL_TUBES = tuple(list(range(29, 59)) + [60])

# Six animals were scanned with only three tissue series; the fourth section was
# too damaged to be worth scanner time.
THREE_SECTION_TUBES = frozenset({30, 33, 34, 42, 53, 54})

# Animals for which a full `_Image_NN_` payload folder is expected to exist.
PAYLOAD_TUBES = (29, 30, 41, 42, 49, 51, 55, 60)

# Preference order when a test needs "some" complete payload dataset.  Tube 49
# is first because it is the four-series animal whose box layout is the
# pipeline's key regression case.
_PAYLOAD_PREFERENCE = (49, 29, 41, 55, 51, 60, 30, 42)


def _read_via_project_resolver() -> Path | None:
    """Ask the pipeline's own `ihc.util.config.load_paths()` where raw data lives.

    Preferred over parsing the YAML here, because that resolver is what the
    pipeline itself uses: it expands `~` and `{root}` templates, so the tests
    and the code under test can never disagree about the data location.
    """
    try:
        from ihc.util.config import load_paths
    except Exception:
        return None
    try:
        paths = load_paths() or {}
    except Exception:
        return None
    if not isinstance(paths, dict):
        return None
    # load_paths() flattens `roots:` to the top level ('raw_root'), but tolerate
    # a nested mapping too in case that changes.
    mappings = [paths]
    nested = paths.get("roots")
    if isinstance(nested, dict):
        mappings.append(nested)
    for mapping in mappings:
        for key in _YAML_ROOT_KEYS:
            value = mapping.get(key)
            if value and isinstance(value, (str, Path)):
                return Path(os.path.expanduser(str(value)))
    return None


def _read_yaml_data_root() -> Path | None:
    """Fallback: read config/paths.yaml directly, without template expansion."""
    if not _PATHS_YAML.is_file():
        return None
    try:
        import yaml
    except ImportError:  # pragma: no cover - pyyaml is a declared dependency
        return None
    try:
        doc = yaml.safe_load(_PATHS_YAML.read_text()) or {}
    except Exception:  # malformed YAML must not be reported as "no config"
        return None
    if not isinstance(doc, dict):
        return None
    candidates = [doc]
    for nested_key in ("roots", "paths"):
        nested = doc.get(nested_key)
        if isinstance(nested, dict):
            candidates.append(nested)
    for mapping in candidates:
        for key in _YAML_ROOT_KEYS:
            value = mapping.get(key)
            # Skip unexpanded templates -- only the project resolver can fill them.
            if isinstance(value, str) and value.strip() and "{" not in value:
                return Path(os.path.expanduser(value.strip()))
    return None


def resolve_data_root() -> tuple[Path | None, str]:
    """Resolve the raw-data root.  Returns (path_or_None, provenance_string)."""
    env = os.environ.get("IHC_DATA_ROOT", "").strip()
    if env:
        return Path(os.path.expanduser(env)), "$IHC_DATA_ROOT"
    from_project = _read_via_project_resolver()
    if from_project is not None:
        return from_project, "ihc.util.config.load_paths() -> roots.raw_root"
    from_yaml = _read_yaml_data_root()
    if from_yaml is not None:
        return from_yaml, str(_PATHS_YAML)
    return DEFAULT_DATA_ROOT, "conftest.DEFAULT_DATA_ROOT fallback"


# --------------------------------------------------------------------------
# dehydration detection (Dropbox online-only placeholders)
# --------------------------------------------------------------------------
def is_materialised(path: Path, *, min_ratio: float = 0.5) -> bool:
    """True when `path` actually occupies disk blocks proportionate to its size.

    A Dropbox "online-only" placeholder reports the real logical size through
    `stat()` but allocates ~0 blocks.  Reading one stalls or errors, so tests
    must detect it up front rather than dying halfway through.
    """
    st = os.stat(path)  # follows symlinks on purpose
    if st.st_size == 0:
        return True
    allocated = getattr(st, "st_blocks", 0) * 512
    return allocated >= st.st_size * min_ratio


def dehydrated_files(root: Path) -> list[Path]:
    """Every regular file under `root` that looks like a cloud placeholder."""
    out = []
    for path in sorted(root.rglob("*")):
        if path.name.startswith("."):
            continue
        if not path.is_file():
            continue
        try:
            if not is_materialised(path):
                out.append(path)
        except OSError:
            out.append(path)
    return out


# --------------------------------------------------------------------------
# markers
# --------------------------------------------------------------------------
def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_data: needs the .vsi index files under the raw-data root",
    )
    config.addinivalue_line(
        "markers",
        "requires_payload: needs a full `_Image_NN_` pixel-payload folder",
    )
    config.addinivalue_line(
        "markers",
        "slow: reads or hashes hundreds of megabytes; deselect with -m 'not slow'",
    )


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
@pytest.fixture(scope="session")
def data_root() -> Path:
    """The raw-data root, or a skip if it is not reachable."""
    root, provenance = resolve_data_root()
    if root is None:
        pytest.skip(
            "No raw-data root configured. Set $IHC_DATA_ROOT or add a "
            f"raw_data: key to {_PATHS_YAML}."
        )
    if not root.is_dir():
        pytest.skip(
            f"Raw-data root {root} (from {provenance}) does not exist. "
            "Set $IHC_DATA_ROOT to override."
        )
    return root


@pytest.fixture(scope="session")
def vsi_paths(data_root: Path) -> dict[int, Path]:
    """{tube_id: path} for every top-level `Image_NN.vsi` index file.

    `Rescan/` is deliberately excluded: it holds re-acquisitions of 51 and 60
    with the same filenames, and folding them in here would silently shadow the
    originals.
    """
    found: dict[int, Path] = {}
    for path in sorted(data_root.glob("Image_*.vsi")):
        stem = path.stem
        digits = stem.split("_", 1)[1] if "_" in stem else ""
        if digits.isdigit():
            found[int(digits)] = path
    if not found:
        pytest.skip(f"No Image_*.vsi index files under {data_root}")
    return found


@pytest.fixture(scope="session")
def payload_datasets(data_root: Path) -> dict[int, tuple[Path, Path]]:
    """{tube_id: (vsi_path, companion_dir)} for animals whose pixels are present.

    Only fully materialised datasets are returned; a dehydrated one is dropped
    here so that downstream tests skip cleanly instead of stalling on a read.
    """
    out: dict[int, tuple[Path, Path]] = {}
    for tube in PAYLOAD_TUBES:
        vsi = data_root / f"Image_{tube}.vsi"
        companion = data_root / f"_Image_{tube}_"
        if not (vsi.is_file() and companion.is_dir()):
            continue
        ets = [p for p in companion.rglob("*.ets")]
        if not ets:
            continue
        if not all(is_materialised(p) for p in [vsi] + ets):
            continue
        out[tube] = (vsi, companion)
    if not out:
        pytest.skip(
            f"No fully materialised `_Image_NN_` payload folder under {data_root}. "
            "Payload folders may be absent or dehydrated (Dropbox online-only)."
        )
    return out


@pytest.fixture(scope="session")
def payload_dataset(payload_datasets) -> tuple[Path, Path]:
    """One complete (vsi_path, companion_dir), preferring tube 49."""
    for tube in _PAYLOAD_PREFERENCE:
        if tube in payload_datasets:
            return payload_datasets[tube]
    return payload_datasets[sorted(payload_datasets)[0]]


@pytest.fixture(scope="session")
def read_meta_cached(vsi_paths):
    """Memoised `read_vsi_meta` so a 31-file sweep parses each index once.

    The import is deliberately lazy: while `ihc.ingest.vsi_meta` is still being
    written, this fixture must not break collection of unrelated test files.
    """
    from ihc.ingest.vsi_meta import read_vsi_meta

    cache: dict[int, object] = {}

    def get(tube: int):
        if tube not in cache:
            path = vsi_paths.get(tube)
            if path is None:
                pytest.skip(f"Image_{tube}.vsi is not present under the data root")
            cache[tube] = read_vsi_meta(path)
        return cache[tube]

    return get
