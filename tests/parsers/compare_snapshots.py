#!/usr/bin/env python3
"""Structural additivity check between two syrupy ``.ambr`` snapshot files.

The migration invariant guarded here is that parser output must never *lose*
archive data. Compared against the pre-migration baseline, the current output
may gain properties and may change how a value is represented, but every piece
of data that existed before must still be present.

The check works in four stages, mirrored by the sections below:

1. **Parse** each ``.ambr`` file into ``{test_name: archive_dict}``.
2. **Flatten** every archive to ``{path: scalar}`` leaves, so two nested
   structures can be compared leaf by leaf.
3. **Diff** the baseline and current leaves per test into three buckets:

   =========  ====================================================  ===========
   bucket     meaning                                               verdict
   =========  ====================================================  ===========
   ``added``   leaf only in current (the newly migrated properties)  fine
   ``changed`` same path, different value (e.g. eV stored as J)       fine
   ``removed`` leaf only in baseline                                  DATA LOSS
   =========  ====================================================  ===========

4. **Report** the buckets and exit non-zero iff anything was removed.

Lists are aligned by index, so a genuine reordering shows up as ``changed``
rather than ``removed``; this is acceptable because the gate only fails on
removals. Deterministic section ordering (from the nomad-lab pre-releases used
in CI) keeps that alignment stable.

Usage: ``compare_snapshots.py <baseline.ambr> <current.ambr>``
"""

from __future__ import annotations

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
def parse_ambr(path: str) -> dict[str, Any]:
    """Parse a syrupy amber file into ``{test_name: value}``.

    The amber format groups each snapshot under a ``# name: <test>`` header,
    followed by an indented body such as ``dict({...})`` / ``list([...])``.
    Those bodies are valid Python once ``dict`` and ``list`` are in scope, so
    they are evaluated in a namespace restricted to those two builtins. The
    files are self-generated CI artifacts, so this is trusted input.
    """
    with open(path) as handle:
        text = handle.read()

    # Split the file into the raw body lines belonging to each snapshot name.
    sections: dict[str, list[str]] = {}
    current: str | None = None
    body: list[str] = []
    for line in text.splitlines():
        if line.startswith('# name:'):
            if current is not None:
                sections[current] = body
            current = line.split(':', 1)[1].strip()
            body = []
        elif line.startswith('#'):
            # '# ---' terminates a snapshot body; other '#' lines are metadata.
            if line.strip() == '# ---' and current is not None:
                sections[current] = body
                current = None
                body = []
        elif current is not None:
            body.append(line)
    if current is not None:
        sections[current] = body

    namespace = {'dict': dict, 'list': list, '__builtins__': {}}
    return {
        name: (eval('\n'.join(lines).strip(), namespace) if lines else None)  # noqa: S307
        for name, lines in sections.items()
    }


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
            return 'DATA LOSS'
        if self.added or self.changed:
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
    baseline = parse_ambr(sys.argv[1])
    current = parse_ambr(sys.argv[2])

    total_removed = 0
    for name in sorted(set(baseline) | set(current)):
        diff = diff_leaves(baseline.get(name), current.get(name))
        total_removed += len(diff.removed)
        report_test(name, baseline.get(name), current.get(name), diff)

    print()
    if total_removed:
        print(f'FAIL: {total_removed} leaves removed relative to baseline (data loss).')
        return 1
    print('PASS: no archive data removed relative to baseline.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
