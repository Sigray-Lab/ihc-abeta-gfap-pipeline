"""Config loading for the whole pipeline.

Owned by the configuration surface. Everything else imports from here rather than
opening YAML or hard-coding a path.

    from ihc.util.config import load_paths, load_config, load_channels, PENDING

`load_paths()` returns a dict of pathlib.Path, with "~" expanded and "{token}"
interpolated against the roots block. `PENDING` is the sentinel string that marks an
undecided value; `pending_items()` finds every one of them.
"""
from pathlib import Path

import yaml

PENDING = "PENDING_PI_DECISION"


def repo_root():
    """The repository root — three parents up from src/ihc/util/config.py."""
    return Path(__file__).resolve().parents[3]


def _read(name):
    with open(repo_root() / "config" / name) as fh:
        return yaml.safe_load(fh)


def load_config():
    """config/config.yaml, unmodified."""
    return _read("config.yaml")


def load_channels():
    """config/channels.yaml, unmodified."""
    return _read("channels.yaml")


def _resolve(value, roots):
    return Path(str(value).format(**roots)).expanduser()


def load_paths():
    """config/paths.yaml with every path resolved to an absolute pathlib.Path.

    Returns a flat dict: top-level roots by name, plus "<section>.<key>" for the
    nested blocks (e.g. "work.qc_dir", "config.slides_csv"). The `raw:` block holds
    file-name patterns rather than paths and is returned verbatim under "raw".
    """
    raw = _read("paths.yaml")
    roots = {"repo_root": str(repo_root())}
    # Ordered so later roots can interpolate earlier ones.
    for key, value in raw["roots"].items():
        if key == "repo_root":
            continue
        roots[key] = str(_resolve(value, roots))

    out = {k: Path(v) for k, v in roots.items()}
    out["raw"] = raw["raw"]
    for section in ("config", "work", "custodian", "results"):
        for key, value in raw.get(section, {}).items():
            out[f"{section}.{key}"] = _resolve(value, roots)
    return out


class PendingDecisionError(RuntimeError):
    """Raised when code tries to use a value the PI has not decided yet.

    This exists because the sentinel string alone is NOT a safety net. It is truthy,
    so `if cfg["some"]["flag"]:` silently takes the yes-branch; it is iterable, so a
    list-shaped pending value iterates as 19 characters; and `random.seed(PENDING)`
    is perfectly valid and yields a deterministic stream. Only the numeric values
    fail loudly on their own.

    So `require()` below is the ONLY sanctioned way to read a config value that
    feeds a calculation. Plain dict access stays available for inventory and
    reporting (`check-config` has to be able to list pending items without raising).
    """


def require(cfg, dotted_key):
    """Read `dotted_key` from `cfg`, refusing to return an undecided value.

    >>> require(load_config(), "gfap_enrichment.bands_um.near_start")

    Raises KeyError if the key is absent, PendingDecisionError if it is still set to
    the PENDING sentinel. Use this everywhere a config value is actually *used*.
    """
    node = cfg
    walked = []
    for part in dotted_key.split("."):
        walked.append(part)
        if not isinstance(node, dict) or part not in node:
            raise KeyError(
                f"config key {dotted_key!r} not found (stopped at {'.'.join(walked)!r})")
        node = node[part]
    if _is_pending(node):
        raise PendingDecisionError(
            f"config key {dotted_key!r} is still {PENDING!r}.\n"
            f"  This is a decision reserved for the PI — see docs/decisions.md.\n"
            f"  The stage that needs it must not run until it is set.")
    return node


def _is_pending(value):
    """True if `value` is the sentinel, or a container holding only sentinels."""
    if isinstance(value, str):
        return value == PENDING
    if isinstance(value, list):
        return bool(value) and all(_is_pending(v) for v in value)
    return False


def pending_items(obj, prefix=""):
    """Every dotted key whose value is the PENDING sentinel, depth-first."""
    found = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            if not prefix and key == "sentinel":
                continue          # the top-level declaration of the token, not a pending value
            found += pending_items(value, f"{prefix}.{key}" if prefix else str(key))
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            found += pending_items(value, f"{prefix}[{i}]")
    elif obj == PENDING:
        found.append(prefix)
    return found


def unapproved_items(obj, prefix=""):
    """Every dotted key whose block carries `pi_approved: false`."""
    found = []
    if isinstance(obj, dict):
        if obj.get("pi_approved") is False:
            found.append(prefix or "<root>")
        for key, value in obj.items():
            found += unapproved_items(value, f"{prefix}.{key}" if prefix else str(key))
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            found += unapproved_items(value, f"{prefix}[{i}]")
    return found
