import json
import os
import re
from importlib import import_module

import click
import numpy as np
import pint


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pint.Quantity):
            return dict(__value=obj.magnitude, __unit=obj.units)
        if isinstance(obj, pint.Unit):
            return str(obj)
        try:
            return super().default(obj)
        except Exception:
            return None


@click.command(help='Generate the output json file for a given parser and mainfile.')
@click.argument('MAINFILE', nargs=1, required=True, type=str)
@click.argument('PARSER_CLS', nargs=1, required=True, type=str)
@click.argument('PARSER_MODULE', nargs=1, required=False, type=str)
def generate_reference_json(mainfile: str, parser_cls: str, parser_module: str = ''):
    if not parser_module:
        match = re.search(r'tests/data/(.+?)/', mainfile)
        if not match:
            raise ValueError(
                'Could not determine parser package from mainfile path. '
                'Please provide the parser module as an argument.'
            )
        package_path = f'nomad_simulation_parsers.parsers.{match.group(1)}.file_parser'
    try:
        package = import_module(package_path)
        parser = getattr(package, parser_cls)()
    except Exception:
        raise ValueError(
            f'Could not load the parser class: {parser_cls} from package: {package}'
        )

    parser.mainfile = mainfile
    parser.parse()

    with open(
        f'{os.path.dirname(mainfile)}/reference_{os.path.basename(mainfile)}.json', 'w'
    ) as f:
        json.dump(parser.to_dict(), f, indent=4, cls=NumpyEncoder)


if __name__ == '__main__':
    generate_reference_json()
