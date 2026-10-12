"""Render parsed blocks for a text parser from the command line."""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Any


def _load_parser_class(import_path: str) -> type[Any]:
    """Load a parser class from a full path or an unqualified class name."""
    module_name, separator, class_name = import_path.rpartition('.')
    if separator:
        try:
            module = importlib.import_module(module_name)
            parser_class = getattr(module, class_name)
        except (ImportError, AttributeError) as error:
            raise ValueError(
                f'Could not import parser {import_path!r}: {error}'
            ) from error
        if not callable(parser_class):
            raise ValueError(f'Parser {import_path!r} is not callable.')
        return parser_class

    package = importlib.import_module('nomad_simulation_parsers.parsers')
    matches = []
    package_path = Path(next(iter(package.__path__)))
    module_names = [package.__name__]
    module_names.extend(
        '.'.join(
            (package.__name__, *path.relative_to(package_path).with_suffix('').parts)
        )
        for path in package_path.rglob('*.py')
        if path.name != '__init__.py'
    )
    for module_name in module_names:
        try:
            parser_class = getattr(importlib.import_module(module_name), import_path)
        except (ImportError, AttributeError):
            continue
        if callable(parser_class):
            matches.append((module_name, parser_class))
    unique_matches = []
    for match in matches:
        if not any(match[1] is existing[1] for existing in unique_matches):
            unique_matches.append(match)
    matches = unique_matches
    if not matches:
        raise ValueError(
            f'Could not find parser {import_path!r} under '
            'nomad_simulation_parsers.parsers.'
        )
    if len(matches) > 1:
        locations = ', '.join(module_name for module_name, _ in matches)
        raise ValueError(f'Parser {import_path!r} is ambiguous; found in {locations}.')
    return matches[0][1]


def _instantiate_parser(parser_class: type[Any], mainfile: Path) -> Any:
    """Instantiate parsers supporting either mainfile= or no arguments."""
    try:
        return parser_class(mainfile=str(mainfile))
    except TypeError as mainfile_error:
        try:
            parser = parser_class()
        except TypeError as constructor_error:
            raise ValueError(
                'Could not instantiate parser with either mainfile= or no '
                f'arguments: {constructor_error}'
            ) from mainfile_error
        if not hasattr(parser, 'mainfile'):
            raise ValueError('Parser instance does not expose a mainfile attribute.')
        parser.mainfile = str(mainfile)
        return parser


def build_argument_parser(
    prog: str = 'nomad-sim-parser view-blocks',
) -> argparse.ArgumentParser:
    """Build the view-blocks command-line parser."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description='Run a TextParser and write an HTML view of its parsed blocks.',
    )
    parser.add_argument(
        'parser',
        metavar='PARSER',
        help='TextParser class name or full import path (package.module.Class).',
    )
    parser.add_argument('mainfile', type=Path, help='Input file passed to the parser.')
    parser.add_argument(
        '--output',
        type=Path,
        help='Destination HTML file (default: <mainfile>.parsed-blocks.html).',
    )
    parser.add_argument(
        '--context-lines',
        type=int,
        default=None,
        help='Show only parsed lines plus this many surrounding lines.',
    )
    parser.add_argument(
        '--key', help='Only show a named quantity or dotted quantity path.'
    )
    parser.add_argument(
        '--leaves-only',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Show only leaf quantities (default: true).',
    )
    return parser


def run(
    arguments: list[str] | None = None,
    *,
    prog: str = 'nomad-sim-parser view-blocks',
) -> int:
    """Run a text parser and write its parsed-block visualization."""
    argument_parser = build_argument_parser(prog)
    args = argument_parser.parse_args(arguments)
    mainfile = args.mainfile.expanduser().resolve()
    if not mainfile.is_file():
        argument_parser.error(f'mainfile is not a file: {mainfile}')
    if args.context_lines is not None and args.context_lines < 0:
        argument_parser.error('--context-lines must be greater than or equal to zero')

    try:
        parser = _instantiate_parser(_load_parser_class(args.parser), mainfile)
        if not hasattr(parser, 'show_visualization'):
            raise ValueError(
                'Parser must provide a show_visualization() method; '
                'use a TextParser subclass.'
            )
        output = args.output or mainfile.with_name(
            f'{mainfile.name}.parsed-blocks.html'
        )
        output = output.expanduser().resolve()
        output = parser.show_visualization(
            context_lines=args.context_lines,
            key=args.key,
            leaves_only=args.leaves_only,
            path=output,
        )
    except (OSError, ValueError) as error:
        print(f'error: {error}', file=sys.stderr)
        return 1

    print(output)
    return 0
