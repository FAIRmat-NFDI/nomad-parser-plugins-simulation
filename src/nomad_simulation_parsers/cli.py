"""Command-line interface for simulation parser utilities."""

from __future__ import annotations

import argparse
import sys

from nomad_simulation_parsers import (
    mapping_report,
    parser_init,
    visualize_parsed_blocks,
)


def _print_main_help() -> None:
    print(
        'usage: nomad-sim-parser '
        '{init,mapping-report,view-blocks} [options]\\n\\n'
        'Commands:\\n'
        '  init             Initialize a parser from metadata YAML.\\n'
        '  mapping-report           Generate the mapping report.\\n'
        '  view-blocks       Run a text parser and write its source view.\\n\\n'
        "Run 'nomad-sim-parser <command> --help' for command options."
    )


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] in {'-h', '--help'}:
        _print_main_help()
        return 0

    command = arguments[0]
    if command == 'mapping-report':
        return mapping_report.run(arguments[1:])
    if command in {'view-blocks', 'visualize-parsed-blocks'}:
        return visualize_parsed_blocks.run(arguments[1:])
    if command == 'init':
        return parser_init.run(arguments[1:])

    parser = argparse.ArgumentParser(prog='nomad-sim-parser')
    parser.error(f'unknown command: {command}')


if __name__ == '__main__':
    raise SystemExit(main())
