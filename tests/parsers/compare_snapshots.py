#!/usr/bin/env python3
"""Structural diff between two parser-output snapshot JSON files.

Compared against the target (baseline) ref, the source (current) ref may only
*add* to the archive: any leaf that existed in the target must still be present
with the same value. Both removals and value modifications are flagged; only
purely additive changes pass.

The check works in four stages, mirrored by the sections below:

1. **Parse** each JSON file into ``{key: archive_dict}``.
2. **Flatten** every archive to ``{path: scalar}`` leaves, so two nested
   structures can be compared leaf by leaf.
3. **Diff** the baseline and current leaves per key into three buckets:

   =========  ====================================================  ===========
   bucket     meaning                                               verdict
   =========  ====================================================  ===========
   ``added``   leaf only in current                                  fine
   ``changed`` same path, different value                            FLAGGED
   ``removed`` leaf only in baseline                                 FLAGGED
   =========  ====================================================  ===========

4. **Report** the buckets and exit non-zero iff anything was removed or changed.

Lists are aligned by index, so a genuine reordering shows up as ``changed``.
Deterministic section ordering (from the nomad-lab pre-releases used in CI)
keeps that alignment stable. Float noise below ``_FLOAT_RTOL`` is ignored.

Usage: ``compare_snapshots.py <target.json> <source.json>``
"""

from __future__ import annotations

import json
import sys
from typing import Any, NamedTuple

# How many example leaves to print per test before truncating.
_MAX_REMOVED_SHOWN = 20
_MAX_CHANGED_SHOWN = 5
# Relative tolerance when comparing two numbers, to ignore float noise.
_FLOAT_RTOL = 1e-9


# --------------------------------------------------------------------------- #
# 1. Parsing
# --------------------------------------------------------------------------- #
def parse_json(path: str) -> dict[str, Any]:
    """Load a snapshot JSON file into ``{key: archive_dict}``."""
    with open(path) as handle:
        return json.load(handle)


# --------------------------------------------------------------------------- #
# 2. Flattening
# --------------------------------------------------------------------------- #
def flatten(
    obj: Any, path: str = '', out: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Flatten a nested archive into ``{path: scalar}`` leaves.

    Dict keys become ``.key`` and list indices become ``[i]`` segments, e.g.
    ``.outputs[0].total_energies[0].value``. Only scalar leaves are stored, so
    two archives can be compared by set operations on their leaf paths.
    """
    if out is None:
        out = {}
    if isinstance(obj, dict):
        for key, value in obj.items():
            flatten(value, f'{path}.{key}', out)
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            flatten(value, f'{path}[{index}]', out)
    elif obj is None and path == '':
        # A wholly-empty archive (archive.data is None) has no leaves, so a
        # None-vs-populated pair reads as purely additive, not a removal.
        pass
    else:
        out[path] = obj
    return out


# --------------------------------------------------------------------------- #
# 3. Diffing
# --------------------------------------------------------------------------- #
class LeafDiff(NamedTuple):
    """Per-test comparison result. ``removed`` being non-empty means data loss."""

    added: list[str]
    removed: list[str]
    changed: list[str]

    @property
    def verdict(self) -> str:
        if self.removed:
            return 'REMOVED'
        if self.changed:
            return 'MODIFIED'
        if self.added:
            return 'additive'
        return 'identical'


def _values_equal(a: Any, b: Any) -> bool:
    """Equality that tolerates floating-point noise between two numbers."""
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= _FLOAT_RTOL * max(1.0, abs(float(a)))
    except (TypeError, ValueError):
        return False


def diff_leaves(baseline: Any, current: Any) -> LeafDiff:
    """Compare one test's baseline and current archives leaf by leaf."""
    base = flatten(baseline)
    cur = flatten(current)
    added = [path for path in cur if path not in base]
    removed = [path for path in base if path not in cur]
    changed = [
        path
        for path in base
        if path in cur and not _values_equal(base[path], cur[path])
    ]
    return LeafDiff(added=added, removed=removed, changed=changed)


# --------------------------------------------------------------------------- #
# 4. Reporting
# --------------------------------------------------------------------------- #
def report_test(name: str, baseline: Any, current: Any, diff: LeafDiff) -> None:
    """Print one line per test plus a few example removed/changed leaves."""
    base = flatten(baseline)
    cur = flatten(current)
    print(
        f'  {diff.verdict:10} {name:34} '
        f'+{len(diff.added)} -{len(diff.removed)} ~{len(diff.changed)}'
    )
    for path in diff.removed[:_MAX_REMOVED_SHOWN]:
        print(f'      removed: {path} = {base[path]!r}')
    for path in diff.changed[:_MAX_CHANGED_SHOWN]:
        print(f'      changed: {path}: {base[path]!r} -> {cur[path]!r}')
    if len(diff.changed) > _MAX_CHANGED_SHOWN:
        print(f'      changed: ... and {len(diff.changed) - _MAX_CHANGED_SHOWN} more')


def main() -> int:
    if len(sys.argv) != 3:
        print(f'usage: {sys.argv[0]} <target.json> <source.json>', file=sys.stderr)
        return 2
    baseline = parse_json(sys.argv[1])
    current = parse_json(sys.argv[2])

    total_removed = total_changed = 0
    for name in sorted(set(baseline) | set(current)):
        diff = diff_leaves(baseline.get(name), current.get(name))
        total_removed += len(diff.removed)
        total_changed += len(diff.changed)
        report_test(name, baseline.get(name), current.get(name), diff)

    print()
    if total_removed or total_changed:
        print(
            f'FAIL: {total_removed} removed, {total_changed} modified relative to the '
            'target. Only additive changes are allowed.'
        )
        return 1
    print('PASS: source only adds to the target (no leaves removed or modified).')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
