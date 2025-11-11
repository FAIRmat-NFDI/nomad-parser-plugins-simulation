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
    
    def get_nmr_text(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        # magnetic shieldings
        data = source.get('ms_list', [])
        from devtools import debug
        debug(data)
        magnetic_shieldings = []
        for atom_data in data:
            values = np.reshape(atom_data[2:], (3, 3))
            FACTOR = 1e-6
            magnetic_shieldings.append(values * FACTOR * ureg("dimensionless"))
        out = dict(magnetic_shieldings=[dict(value=m) for m in magnetic_shieldings ])

        # magnetic_susceptibilities
        chi_bare_pGv = source.get("chi_bare_pGv", [])
        chi_bare_vGv = source.get("chi_bare_vGv", [])

        sus = (chi_bare_pGv + chi_bare_vGv) / 2
        out["magnetic_susceptibilities"] = dict(
            value=sus,
            value_vgv_approx = chi_bare_vGv,
            value_pgv_approx = chi_bare_pGv
            )
        debug(out)
        
        
        return [out]



class GIPAWMainfileXMLParser(MainfileXMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER
    
    def get_nmr_xml(self, atom: dict[str, Any], **kwargs) -> Any:
        value = np.reshape(atom["__value"], (3, 3))
        FACTOR = 1e-6
        return value * FACTOR * ureg("dimensionless")



class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    _text_parser = GIPAWMainfileTextParser(text_parser=GIPAWFileParser())
    _xml_parser = GIPAWMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
