"""Load config/config.yaml and resolve output paths against the repository root."""
import os
from functools import lru_cache
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "config" / "config.yaml"


@lru_cache(maxsize=None)
def load_config(path=None):
    """Read the YAML config. `GXE_CONFIG` overrides the default location."""
    path = Path(path or os.environ.get("GXE_CONFIG", DEFAULT_CONFIG))
    with open(path) as fh:
        return yaml.safe_load(fh)


def cohort_names(cfg=None):
    cfg = cfg or load_config()
    return list(cfg["cohorts"])


def _resolve(p):
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


def intermediate_path(*parts, cfg=None):
    """Path under the intermediate data dir; parent directories are created."""
    cfg = cfg or load_config()
    out = _resolve(cfg["outputs"]["intermediate_dir"]).joinpath(*parts)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def assembled_path(*parts, cfg=None):
    cfg = cfg or load_config()
    out = _resolve(cfg["outputs"]["assembled_dir"]).joinpath(*parts)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out
