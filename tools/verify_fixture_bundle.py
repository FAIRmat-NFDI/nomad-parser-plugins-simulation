"""Verify that nightly-test fixtures match their committed SHA-256 manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

SHA256_LENGTH = 64


def parse_manifest(path: Path) -> list[tuple[str, Path]]:
    entries = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            digest, filename = line.split(maxsplit=1)
        except ValueError as error:
            msg = f'{path}:{line_number}: expected "SHA256 relative/path"'
            raise ValueError(msg) from error
        if len(digest) != SHA256_LENGTH or any(
            char not in '0123456789abcdef' for char in digest
        ):
            msg = f'{path}:{line_number}: invalid SHA-256 digest'
            raise ValueError(msg)
        entries.append((digest, Path(filename)))
    if not entries:
        msg = f'{path}: fixture manifest must not be empty'
        raise ValueError(msg)
    return entries


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('manifest', type=Path)
    args = argument_parser.parse_args()

    for expected, relative_path in parse_manifest(args.manifest):
        if relative_path.is_absolute() or '..' in relative_path.parts:
            raise ValueError(f'{args.manifest}: unsafe fixture path {relative_path}')
        if not relative_path.is_file():
            raise FileNotFoundError(f'missing required fixture: {relative_path}')
        actual = sha256(relative_path)
        if actual != expected:
            msg = (
                f'checksum mismatch for {relative_path}: expected {expected}, '
                f'got {actual}'
            )
            raise ValueError(msg)


if __name__ == '__main__':
    main()
