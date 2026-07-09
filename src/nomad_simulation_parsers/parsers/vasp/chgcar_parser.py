import os
import re

import numpy as np
from nomad.parsing.file_parser.file_parser import FileParser
from nomad.parsing.file_parser.mapping_parser import MappingParser
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.utils.general import search_files

LOGGER = get_logger(__name__)


class CHGCARFileParser(FileParser):
    def parse(self, key=None):
        if self._results is None:
            self._results = {}

        values = []
        with self.open_mainfile_obj() as f:
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
                    values.append(
                        np.reshape(np.array(charge_density, np.float64), grid)
                    )
                    grid = []
                    n_points = 0
                    charge_density = []
        self._results['values'] = values


class CHGCARParser(MappingParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs) -> dict:
        if self.data_object:
            self.data_object.parse()
            return self.data_object._results
        return {}

    def load_file(self) -> FileParser | None:
        chgcar_files = search_files('*CHGCAR*', os.path.dirname(self.filepath))
        if not chgcar_files:
            return None

        if len(chgcar_files) > 1:
            self.logger.warning('Found more than one CHGCAR file, parsing only first.')
        chgcar_parser = CHGCARFileParser()
        chgcar_parser.mainfile = chgcar_files[0]
        return chgcar_parser

    def from_dict(self, dct: dict):
        pass
