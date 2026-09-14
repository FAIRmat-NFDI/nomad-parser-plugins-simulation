import re
from typing import Any

from nomad_file_parser.text_parser import Quantity, TextParser

from .common import to_float


def str_to_input_parameters(block: str) -> dict[str, Any]:
    re_array = re.compile(r'\s*([\w\-]+)\[[\d ]+\]\s*=\s*\{*(.+)')
    re_scalar = re.compile(r'\s*([\w\-]+)\s*[=:]\s*(.+)')
    parameters = dict()
    val = [line.strip() for line in block.splitlines()]
    for val_n in val:
        if match := re_scalar.match(val_n):
            parameters[match.group(1)] = match.group(2)
        elif match := re_array.match(val_n):
            parameters.setdefault(match.group(1), [])
            value = [to_float(v) for v in match.group(2).rstrip('}').split(',')]
            parameters[match.group(1)].append(value[0] if len(value) == 1 else value)
    return parameters


class GromacsMdpParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = [
            Quantity(
                'input_parameters',
                r'([\s\S]+)',
                str_operation=str_to_input_parameters,
            ),
        ]
