"""Command-line interface for simulation parser utilities."""

from __future__ import annotations

import argparse
import sys

from nomad_simulation_parsers import mapping_report
from nomad_simulation_parsers import parser_init


def _print_main_help() -> None:
    print(
        'usage: nomad-sim-parser {init,mapping-report} [options]\\n\\n'
        'Commands:\\n'
        '  init             Initialize a parser from metadata YAML.\\n'
        '  mapping-report   Generate the mapping report.\\n\\n'
        "Run 'nomad-sim-parser <command> --help' for command options."
    )


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] in {'-h', '--help'}:
        _print_main_help()
        return 0

    command = arguments[0]
    if command not in {'init', 'mapping-report'}:
        parser = argparse.ArgumentParser(prog='nomad-sim-parser')
        parser.error(f'unknown command: {command}')

    if command == 'mapping-report':
        return mapping_report.run(arguments[1:])

    return parser_init.run(arguments[1:])


if __name__ == '__main__':
    raise SystemExit(main())
