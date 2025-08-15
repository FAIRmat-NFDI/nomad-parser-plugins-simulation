import os
import re

import numpy as np
from nomad.parsing.file_parser.mapping_parser import MappingParser
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.utils.general import search_files

LOGGER = get_logger(__name__)


class CHGCARParser(MappingParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs):
        dct = dict(values=[])
        with open(self.filepath) as f:
            grid = None
            n_points = 0
            charge_density = []
            re_grid = re.compile(r' *\d+ +\d+ +\d+\s+')
            N = -1
            for line in f:
                N += 1
                if not line.strip():
                    grid = []
                if grid is None:
                    continue

                match = re_grid.match(line)
                if match and n_points == 0:
                    grid = [int(i) for i in line.strip().split()]
                    n_points = grid[0] * grid[1] * grid[2]
                elif len(charge_density) < n_points:
                    charge_density.extend([float(v) for v in line.strip().split()])
                if charge_density and len(charge_density) == n_points:
                    dct['values'].append(
                        np.reshape(np.array(charge_density, np.float64), grid)
                    )
                    grid = []
                    n_points = 0
                    charge_density = []
        return dct

    def load_file(self) -> dict:
        return {}

    def from_dict(self, dct: dict):
        pass


def parse_chgcar(chgcar_file: str, archive_parser: MappingParser) -> None:
    if not archive_parser.data_object.m_root().m_context:
        return

    chgcar_files = search_files(
        os.path.basename(chgcar_file), os.path.dirname(chgcar_file)
    )
    if not chgcar_files:
        return

    chgcar_parser = CHGCARParser()
    if len(chgcar_files) > 1:
        chgcar_parser.logger.warning(
            'Found more than one CHGCAR file, parsing only first.'
        )
    chgcar_parser.filepath = chgcar_files[0]
    chgcar_parser.convert(archive_parser)
