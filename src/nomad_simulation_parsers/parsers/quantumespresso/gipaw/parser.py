from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad.units import ureg

from nomad_simulation_parsers.schema_packages.quantumespresso import gipaw

from ..parser import MainfileTextParser, MainfileXMLParser, QuantumEspressoArchiveWriter
from .file_parser import GIPAWFileParser

LOGGER = get_logger(__name__)


class GIPAWMainfileTextParser(MainfileTextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER
    
    def get_magnetic_shieldings(self, source: dict[str, Any]) -> list[dict[str, Any]]:

        data = source.get('ms_list', [])

        magnetic_shieldings = []
        for atom_data in data:
            values = np.reshape(atom_data[2:], (3, 3))
            FACTOR = 1e-6
            magnetic_shieldings.append(values * FACTOR * ureg("dimensionless"))
        return [dict(magnetic_shieldings=[dict(value=m) for m in magnetic_shieldings ])]

        



class GIPAWMainfileXMLParser(MainfileXMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    _text_parser = GIPAWMainfileTextParser(text_parser=GIPAWFileParser())
    _xml_parser = GIPAWMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
