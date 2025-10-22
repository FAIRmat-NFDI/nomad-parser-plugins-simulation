import re
from typing import Any

from nomad.parsing.file_parser.text_parser import Quantity, TextParser

from .common import RE_N, to_float


def str_to_header(block: str) -> dict[str, Any]:
    n_val = 2
    val = [v.split(':', 1) for v in block.strip().splitlines()]
    return {v[0].strip(): v[1].strip() for v in val if len(v) == n_val}


def str_to_input_parameters(block: str) -> dict[str, Any]:
    re_section = re.compile(r'^\s*([\w\-]+):\s*$')
    re_subsection = re.compile(r'^\s*([\w\-]+\s[\d]+):\s*$')
    re_scalar = re.compile(r'\s*([\w\-]+)\s*[=:]\s*(.+)')
    re_array = re.compile(r'\s*([\w\-]+)\[[\d ]+\]\s*=\s*\{*(.+)')
    re_shorthand_array = re.compile(
        r'\s*([\w\-]+)\[\d+,\.\.\.,\d+\]\s*=\s*\{(\d+),\.\.\.,(\d+)\}'
    )

    parameters = dict()
    stack = [parameters]  # Stack to track the current context
    indent_levels = []  # To track the indentation levels

    for line in block.strip().splitlines():
        val_n = line.rstrip()  # Remove trailing spaces
        if not val_n:
            continue

        # Calculate the indentation level
        current_indent = len(val_n) - len(val_n.lstrip())

        # Handle end of section based on indentation
        while indent_levels and current_indent <= indent_levels[-1]:
            stack.pop()
            indent_levels.pop()

        # Check for section or subsection
        if match := (re_section.match(val_n) or re_subsection.match(val_n)):
            key = match.group(1)
            stack[-1][key] = {}
            stack.append(stack[-1][key])
            indent_levels.append(current_indent)
        # Check for scalar
        elif match := re_scalar.match(val_n):
            key = match.group(1)
            value = match.group(2)
            if value.lower() in ['true', 'false']:
                value = value.lower() == 'true'
            elif value.replace('.', '', 1).isdigit():
                value = float(value) if '.' in value else int(value)
            stack[-1][key] = value
        # Check for shorthand array
        elif match := re_shorthand_array.match(val_n):
            array_key = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            stack[-1][array_key] = list(range(start, end + 1))
        # Check for array
        elif match := re_array.match(val_n):
            array_key = match.group(1)
            value = [float(v) for v in match.group(2).rstrip('}').split(',')]
            stack[-1].setdefault(array_key, [])
            stack[-1][array_key].append(value[0] if len(value) == 1 else value)
    return parameters


def str_to_energies(block: str) -> dict[str, float]:
    thermo_common = [
        r'Total Energy',
        r'Potential',
        r'Kinetic En.',
        r'Temperature',
        r'Pressure \(bar\)',
        r'LJ \(SR\)',
        r'Coulomb \(SR\)',
        r'Proper Dih.',
    ]
    n_chars_val = re.search(rf'( +{"| +".join(thermo_common)})', block)
    n_chars_val = len(n_chars_val.group(1)) if n_chars_val is not None else None
    if n_chars_val is None:
        n_chars_val = 15
    energies = {}
    rows = [v for v in block.splitlines() if v]
    for n in range(0, len(rows), 2):
        pointer = 0
        while pointer < len(rows[n]):
            key = rows[n][pointer : pointer + n_chars_val].strip()
            value = rows[n + 1][pointer : pointer + n_chars_val]
            float_value = to_float(value)
            if float_value is not None:
                energies[key] = to_float(value)
            pointer += n_chars_val
    return energies


def str_to_step_info(block: str) -> dict[str, float]:
    val = block.strip().splitlines()
    keys = val[0].split()
    values = [to_float(v) for v in val[1].split()]
    return {key: values[n] for n, key in enumerate(keys) if values[n] is not None}


class GromacsLogParser(TextParser):
    def init_quantities(self):
        thermo_quantities = [
            Quantity(
                'energies',
                r'Energies \(kJ/mol\).*\n(\s*[\s\S]+?)(?:\n.*step.* load imb.*|\n\n)',
                str_operation=str_to_energies,
                convert=False,
            ),
            Quantity(
                'step_info',
                rf'{RE_N}\s*(Step.+\n[\d\.\- ]+)',
                str_operation=str_to_step_info,
                convert=False,
            ),
        ]

        self._quantities = [
            Quantity('time_start', r'Log file opened on (.+)', flatten=False),
            Quantity(
                'host_info',
                r'Host:\s*(\S+)\s*pid:\s*(\d+)\s*rank ID:\s*(\d+)\s*number of ranks:\s*(\d*)',
            ),
            Quantity(
                'module_version', r'GROMACS:\s*(.+?),\s*VERSION\s*(\S+)', flatten=False
            ),
            Quantity('execution_path', r'Executable:\s*(.+)'),
            Quantity('working_path', r'Data prefix:\s*(.+)'),
            # TODO cannot understand treatment of the command line in the old parser
            Quantity(
                'header',
                r'(?:GROMACS|Gromacs) (20[\s\S]+?)\n\n',
                str_operation=str_to_header,
            ),
            Quantity(
                'header',
                r'(?:GROMACS|Gromacs) (version:[\s\S]+?)\n\n',
                str_operation=str_to_header,
            ),
            Quantity(
                'input_parameters',
                r'Input Parameters:\s*\n([\s\S]+?)\n\n',
                str_operation=str_to_input_parameters,
            ),
            Quantity('maximum_force', r'Norm of force\s*([\s\S]+?)\n\n', flatten=False),
            Quantity(
                'step',
                r'(Step\s*Time[\s\S]+?Energies[\s\S]+?\n\n)',
                repeats=True,
                sub_parser=TextParser(quantities=thermo_quantities),
            ),
            Quantity(
                'averages',
                r'A V E R A G E S  ====>([\s\S]+?\n\n\n)',
                sub_parser=TextParser(quantities=thermo_quantities),
            ),
            Quantity('time_end', r'Finished \S+ on rank \d+ (.+)', flatten=False),
        ]
