#!/usr/bin/env python3
"""Structural diff between two parser-output snapshot JSON files.

Compared against the target (baseline) ref, the source (current) ref may only
*add* to the archive: any leaf that existed in the target must still be present
with the same value. This yields three tiers of signal per parser:

===========  =====================================  ==========================
verdict      meaning                                outcome
===========  =====================================  ==========================
``identical`` no leaf added, removed, or changed     pass, silent
``additive``  only additions                         pass, surfaced as a
                                                     non-failing GitHub warning
``MODIFIED`` / a leaf changed value or was removed   hard fail (exit 1),
``REMOVED``                                          printing ``old -> new``
===========  =====================================  ==========================

So additions stay visible without turning the check red, while any removal or
value change blocks. This policy is intentional: a deliberate move or value
change is meant to fail here, and is reviewed by running the check locally.

The check works in stages, mirrored by the sections below:

1. **Parse** each JSON file into its ``provenance`` block and ``{key: archive}``
   snapshots. If the two snapshot sets are byte-identical the diff is skipped
   entirely (a fast hash short-circuit). Provenance is reported for context but
   never diffed -- two refs may legitimately differ there.
2. **Flatten** every archive to ``{path: scalar}`` leaves, so two nested
   structures can be compared leaf by leaf.
3. **Diff** the baseline and current leaves per key into added / changed /
   removed buckets.
4. **Report** the buckets (and emit warning annotations for additive parsers
   under GitHub Actions) and exit non-zero iff anything was removed or changed.

Lists are aligned by index, so a genuine reordering shows up as ``changed``.
Deterministic section ordering (from the nomad-lab pre-releases used in CI)
keeps that alignment stable. Float noise below ``_FLOAT_RTOL`` is ignored.

Usage: ``compare_snapshots.py <target.json> <source.json>``
"""

from __future__ import annotations

import json
import os
import sys
from typing import Any, NamedTuple

# How many example leaves to print per test before truncating.
_MAX_ADDED_SHOWN = 5
_MAX_REMOVED_SHOWN = 20
_MAX_CHANGED_SHOWN = 5
# Relative tolerance when comparing two numbers, to ignore float noise.
_FLOAT_RTOL = 1e-9
# True when running inside a GitHub Actions job, where ``::warning::`` lines
# render as non-failing annotations in the checks UI.
_ON_GITHUB = os.environ.get('GITHUB_ACTIONS') == 'true'


# --------------------------------------------------------------------------- #
# 1. Parsing
# --------------------------------------------------------------------------- #
def parse_json(path: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load a snapshot JSON file into ``(snapshots, provenance)``.

    Accepts the current ``{"provenance": ..., "snapshots": ...}`` layout and,
    for robustness, a bare ``{key: archive}`` mapping (provenance then empty).
    """
    with open(path) as handle:
        document = json.load(handle)
    if isinstance(document, dict) and 'snapshots' in document:
        return document['snapshots'], document.get('provenance', {})
    return document, {}


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
    # Sort each bucket so the report ordering is stable across runs regardless of
    # dict iteration order; this does not affect pass/fail semantics.
    added = sorted(path for path in cur if path not in base)
    removed = sorted(path for path in base if path not in cur)
    changed = sorted(
        path
        for path in base
        if path in cur and not _values_equal(base[path], cur[path])
    )
    return LeafDiff(added=added, removed=removed, changed=changed)


# --------------------------------------------------------------------------- #
# 4. Reporting
# --------------------------------------------------------------------------- #
def report_test(name: str, baseline: Any, current: Any, diff: LeafDiff) -> None:
    """Print one line per test plus a few example added/removed/changed leaves.

    Additive parsers (additions only, nothing removed or changed) also emit a
    non-failing ``::warning::`` annotation under GitHub Actions, so growth is
    visible in the checks UI without turning the job red.
    """
    base = flatten(baseline)
    cur = flatten(current)
    print(
        f'  {diff.verdict:10} {name:34} '
        f'+{len(diff.added)} -{len(diff.removed)} ~{len(diff.changed)}'
    )
    for path in diff.added[:_MAX_ADDED_SHOWN]:
        print(f'      added:   {path} = {cur[path]!r}')
    if len(diff.added) > _MAX_ADDED_SHOWN:
        print(f'      added:   ... and {len(diff.added) - _MAX_ADDED_SHOWN} more')
    for path in diff.removed[:_MAX_REMOVED_SHOWN]:
        print(f'      removed: {path} = {base[path]!r}')
    for path in diff.changed[:_MAX_CHANGED_SHOWN]:
        print(f'      changed: {path}: {base[path]!r} -> {cur[path]!r}')
    if len(diff.changed) > _MAX_CHANGED_SHOWN:
        print(f'      changed: ... and {len(diff.changed) - _MAX_CHANGED_SHOWN} more')

    if _ON_GITHUB and diff.verdict == 'additive':
        sample = ', '.join(diff.added[:_MAX_ADDED_SHOWN])
        print(
            f'::warning title=additive-output-changes::{name}: '
            f'{len(diff.added)} leaf(s) added (e.g. {sample})'
        )


def _print_provenance(label: str, prov: dict[str, Any]) -> None:
    if not prov:
        return
    fixtures = prov.get('fixtures') or {}
    print(
        f'  {label}: python {prov.get("python")}, nomad-lab {prov.get("nomad_lab")}, '
        f'uv.lock {str(prov.get("uv_lock_sha256"))[:12]}, {len(fixtures)} fixture(s)'
    )


def main() -> int:
    if len(sys.argv) != 3:
        print(f'usage: {sys.argv[0]} <target.json> <source.json>', file=sys.stderr)
        return 2
    baseline, baseline_prov = parse_json(sys.argv[1])
    current, current_prov = parse_json(sys.argv[2])

    print('provenance:')
    _print_provenance('target', baseline_prov)
    _print_provenance('source', current_prov)
    print()

    # Fast hash short-circuit: if the two snapshot sets are byte-identical there
    # is nothing to diff, so skip straight to a pass.
    if json.dumps(baseline, sort_keys=True) == json.dumps(current, sort_keys=True):
        print('PASS: source snapshots are identical to the target (hash match).')
        return 0

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
