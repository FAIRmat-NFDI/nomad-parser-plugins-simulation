#!/usr/bin/env python3
"""Structural additivity check between two syrupy ``.ambr`` snapshot files.

The migration invariant is that parser output must never *lose* archive data:
every leaf present in the pre-migration baseline must still be present in the
current output (values may change representation, dicts may gain keys, lists may
grow). This script flattens both snapshots to ``path -> value`` leaves, aligns
lists by index, and reports per test:

* ``added``   - leaves only in current (expected: the migrated properties)
* ``changed`` - same path, different value (e.g. eV->J unit storage)
* ``removed`` - leaves only in baseline == DATA LOSS

Exit code is non-zero iff any leaf was removed.

Usage: ``compare_snapshots.py <baseline.ambr> <current.ambr>``
"""

from __future__ import annotations

import sys


def parse_ambr(path: str) -> dict:
    """Parse a syrupy amber file into ``{test_name: value}``.

    Bodies use ``dict({...})`` / ``list([...])`` which are evaluated in a
    restricted namespace (the files are self-generated, trusted CI artifacts).
    """
    with open(path) as handle:
        text = handle.read()
    sections: dict[str, list[str]] = {}
    current: str | None = None
    buffer: list[str] = []
    for line in text.splitlines():
        if line.startswith('# name:'):
            if current is not None:
                sections[current] = buffer
            current = line.split(':', 1)[1].strip()
            buffer = []
        elif line.startswith('#'):
            if line.strip() == '# ---' and current is not None:
                sections[current] = buffer
                current = None
                buffer = []
        elif current is not None:
            buffer.append(line)
    if current is not None:
        sections[current] = buffer

    namespace = {'dict': dict, 'list': list, '__builtins__': {}}
    result = {}
    for name, lines in sections.items():
        body = '\n'.join(lines).strip()
        result[name] = eval(body, namespace) if body else None  # noqa: S307
    return result


def flatten(obj, path: str = '', out: dict | None = None) -> dict:
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


def _close(a, b) -> bool:
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(a)))
    except (TypeError, ValueError):
        return False


def main() -> int:
    baseline = parse_ambr(sys.argv[1])
    current = parse_ambr(sys.argv[2])

    total_removed = 0
    for test in sorted(set(baseline) | set(current)):
        base_leaves = flatten(baseline.get(test))
        cur_leaves = flatten(current.get(test))
        added = [p for p in cur_leaves if p not in base_leaves]
        removed = [p for p in base_leaves if p not in cur_leaves]
        changed = [
            p
            for p in base_leaves
            if p in cur_leaves and not _close(base_leaves[p], cur_leaves[p])
        ]
        total_removed += len(removed)
        status = 'DATA LOSS' if removed else ('additive' if added or changed else 'identical')
        print(
            f'  {status:10} {test:34} '
            f'+{len(added)} -{len(removed)} ~{len(changed)}'
        )
        for path in removed[:20]:
            print(f'      removed: {path} = {base_leaves[path]!r}')
        for path in changed[:5]:
            print(
                f'      changed: {path}: '
                f'{base_leaves[path]!r} -> {cur_leaves[path]!r}'
            )
        if len(changed) > 5:
            print(f'      changed: ... and {len(changed) - 5} more')

    print()
    if total_removed:
        print(f'FAIL: {total_removed} leaves removed relative to baseline (data loss).')
        return 1
    print('PASS: no archive data removed relative to baseline.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
