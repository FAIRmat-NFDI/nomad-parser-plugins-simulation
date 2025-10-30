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
        from devtools import debug
        debug("sto facendo qualcosa di nuovo")
        data = source.get('ms_list', [])

        tensors = np.array(
            [np.array(row[2:], dtype=np.float64).reshape(3, 3) for row in data],
            dtype=np.float64,
        )
        FACTOR = 1e-6
        magnetic_shieldings = (tensors * FACTOR) * ureg("dimensionless")
        result = [dict(magnetic_shieldings=magnetic_shieldings)]
        return result
        



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
