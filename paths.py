"""
Load and resolve paths.yaml.

Used by every stage that needs configured paths. Resolves `${var}` references
recursively against other top-level keys; values written as absolute paths are
returned unchanged, so users can point any single key outside `work_dir`.

Lookup order for the yaml file:
  1. explicit `yaml_path` arg to `load_paths`
  2. `IMPACT_DL_PATHS` env var
  3. `<repo_root>/paths.yaml` (alongside this module)
"""

import os
import re
from pathlib import Path

import yaml

_VAR_PATTERN = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}")
_REPO_ROOT = Path(__file__).resolve().parent


def _resolve(raw, key, seen):
    if key in seen:
        raise ValueError("Cyclic reference in paths.yaml at '{}'".format(key))
    val = raw[key]
    if not isinstance(val, str):
        return val

    def sub(m):
        ref = m.group(1)
        if ref not in raw:
            raise KeyError("paths.yaml refers to undefined key '${{{}}}'".format(ref))
        resolved = _resolve(raw, ref, seen | {key})
        return str(resolved)

    return _VAR_PATTERN.sub(sub, val)


def load_paths(yaml_path=None, **overrides):
    """
    Load paths.yaml, resolve `${var}` interpolations, and apply CLI overrides.

    Returns a plain dict (so it works as a Snakemake config).
    """
    if yaml_path is None:
        yaml_path = os.environ.get("IMPACT_DL_PATHS") or (_REPO_ROOT / "paths.yaml")
    yaml_path = Path(yaml_path).expanduser()

    if not yaml_path.exists():
        raise FileNotFoundError(
            f"paths.yaml not found at {yaml_path}. "
            f"Copy paths.yaml.example to paths.yaml and edit it."
        )

    with open(yaml_path) as f:
        raw = yaml.safe_load(f)

    resolved = {k: _resolve(raw, k, set()) for k in raw}

    for k, v in overrides.items():
        if v is not None:
            resolved[k] = v

    return resolved


if __name__ == "__main__":
    import argparse
    import json
    import sys

    parser = argparse.ArgumentParser(description="Resolve and print paths.yaml")
    parser.add_argument("--paths_yaml", default=None)
    parser.add_argument("--get", default=None,
                        help="Print a single resolved value (for bash scripts).")
    args = parser.parse_args()
    paths = load_paths(args.paths_yaml)
    if args.get:
        if args.get not in paths:
            print(f"ERROR: key '{args.get}' not in paths.yaml", file=sys.stderr)
            sys.exit(1)
        print(paths[args.get])
    else:
        print(json.dumps(paths, indent=2, default=str))
