#!/usr/bin/env python3
"""Reject a workflow GitHub would refuse to parse.

`release.yml` runs on a tag and on nothing else, so the first time anybody
learns it is malformed is the moment of a release. That happened: an `exclude`
entry named `target`, a key that existed only inside `include`, and GitHub
rejected the whole file before a single job started. No test saw it, no review
caught it, and the other CI jobs cannot: they never read this file.

Two checks, both of rules GitHub enforces and neither of which actionlint
covers (verified against 1.7.7, which exits 0 on the file GitHub refused).

1. Every workflow is valid YAML.
2. Every key named in a matrix `exclude` is a matrix *axis*. GitHub builds the
   cross product from the axes alone, applies `exclude` to that, and only then
   merges `include`. So a key contributed by `include` is not one `exclude` can
   name, and naming it is a parse error rather than a no-op.

Run with no arguments from the repository root.
"""

from __future__ import annotations

import pathlib
import sys

import yaml

# Keys that live beside the axes in `matrix:` and are not axes themselves.
NOT_AN_AXIS = {"include", "exclude"}


def problems(path: pathlib.Path) -> list[str]:
    try:
        doc = yaml.safe_load(path.read_text())
    except yaml.YAMLError as e:
        return [f"{path}: not valid YAML: {e}"]
    if not isinstance(doc, dict):
        return [f"{path}: top level is not a mapping"]

    found = []
    for job_name, job in (doc.get("jobs") or {}).items():
        if not isinstance(job, dict):
            continue
        matrix = (job.get("strategy") or {}).get("matrix")
        if not isinstance(matrix, dict):
            continue
        axes = set(matrix) - NOT_AN_AXIS
        for entry in matrix.get("exclude") or []:
            if not isinstance(entry, dict):
                continue
            for key in entry:
                if key not in axes:
                    found.append(
                        f"{path}: job '{job_name}': matrix exclude names "
                        f"'{key}', which is not a matrix axis "
                        f"({', '.join(sorted(axes)) or 'none'}). GitHub "
                        f"refuses the whole workflow for this. Make '{key}' "
                        f"an axis, and attach the rest with include."
                    )
    return found


def main() -> int:
    root = pathlib.Path(__file__).resolve().parent.parent / ".github" / "workflows"
    files = sorted(p for p in root.iterdir() if p.suffix in (".yml", ".yaml"))
    if not files:
        print(f"check_workflows: no workflows under {root}", file=sys.stderr)
        return 1
    found = [msg for f in files for msg in problems(f)]
    for msg in found:
        print(f"error: {msg}", file=sys.stderr)
    if found:
        return 1
    print(f"check_workflows: ok, {len(files)} workflow(s) parse")
    return 0


if __name__ == "__main__":
    sys.exit(main())
