from typing import Any

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import pwscf

from ..parser import MainfileParser
from .file_parser import PWSCFFileParser

LOGGER = get_logger(__name__)


class PWSCFMainfileParser(MainfileParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_configurations(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        methods = {
            'self_consistent': 'single_point',
            'bandstructure': 'single_point',
            'bfgs_geometry_optimization': 'geometry_optimization',
            'molecular_dynamics': 'molecular_dynamics',
            'langevin_dynamics': 'langevin_dynamics',
            'damped_dynamics': 'geometry_optimization',
            'vcs_wentzcovitch_damped_minimization': 'geometry_optimization',
        }

        configurations = []
        for key in methods:
            config = source.get(key)
            if config is None:
                continue
            configurations.append(config.get('self_consistent', config))
        return configurations


class PWSCFArchiveWriter(QuantumEspressoArchiveWriter):
    schema = pwscf
    mainfile_parser = PWSCFMainfileParser(text_parser=PWSCFFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = 'out'
        super().parse_program(archive, index)
