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
        from devtools import debug
        debug("sono get_nmr_text")
        out = {}
        # magnetic shieldings
        ms_list = source.get('ms_list', None)
        if ms_list is not None:
            magnetic_shieldings = []
            for atom_data in ms_list:
                values = np.reshape(atom_data[2:], (3, 3))
                FACTOR = 1e-6
                magnetic_shieldings.append(values * FACTOR * ureg("dimensionless"))
            out["magnetic_shieldings"] = [dict(value=m) for m in magnetic_shieldings ]

        # magnetic_susceptibilities
        chi_bare_pGv = source.get("chi_bare_pGv", None)
        chi_bare_vGv = source.get("chi_bare_vGv", None)

        if chi_bare_pGv is not None and chi_bare_pGv is not None:
            sus = (chi_bare_pGv + chi_bare_vGv) / 2
            out["magnetic_susceptibilities"] = dict(
                value=sus,
                value_vgv_approx = chi_bare_vGv,
                value_pgv_approx = chi_bare_pGv
                )

        # electric field gradient
        efg = source.get('efg', None)
        if efg is not None:
            electric_field_gradients = []
            for i, atom_data in enumerate(efg):
                values = np.reshape(atom_data[2:], (3, 3))
                electric_field_gradients.append(values)
            out['electric_field_gradients'] = [dict(value=e) for e in electric_field_gradients]
        
        return [out]



class GIPAWMainfileXMLParser(MainfileXMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER
    
    def get_magnetic_shieldings(self, atom: dict[str, Any], **kwargs) -> Any:
        value = np.reshape(atom["__value"], (3, 3))
        FACTOR = 1e-6
        return value * FACTOR * ureg("dimensionless")
    

    def get_magnetic_susceptibilities(self, source: dict[str, Any], **kwargs) -> Any:
        if kwargs["name"] != "value":
            value = source.get(kwargs["name"], None)
            return np.reshape(value.get("__value", None), (3, 3))

        value_vgv = source.get("susceptibility_low", None)
        vgv = np.reshape(value_vgv.get("__value", None), (3, 3))
        value_pgv = source.get("susceptibility_high", None)
        pgv = np.reshape(value_pgv.get("__value", None), (3, 3))
        sus = (vgv + pgv) / 2
        return sus
        



class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    _text_parser = GIPAWMainfileTextParser(text_parser=GIPAWFileParser())
    _xml_parser = GIPAWMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
