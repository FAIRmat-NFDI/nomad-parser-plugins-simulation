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

        # electric_field_gradient
        data = source.get('efg', None)
        if data is not None:
            electric_field_gradients = []
            for atom_data in data:
                values = np.reshape(atom_data[2:], (3, 3))
                electric_field_gradients.append(values)
            out['electric_field_gradients'] = [dict(value=item) for item in electric_field_gradients]

        # hyperfine_dipolar
        data = source.get('hyperfine_dipolar', None)
        if data is not None:
            hyperfine_dipolar = []
            for atom_data in data:
                values = np.reshape(atom_data[2:], (3, 3))
                hyperfine_dipolar.append(values)
            out['hyperfine_dipolar'] = [dict(value=item) for item in hyperfine_dipolar]

        # hyperfine_fermi_contact
        data = source.get('hyperfine_fermi_contact', None)
        if data is not None:
            hyperfine_fermi_contact = []
            for atom_data in data:
                values = atom_data[-1]
                hyperfine_fermi_contact.append(values)
            out['hyperfine_fermi_contact'] = [dict(value=item) for item in hyperfine_fermi_contact]

        # delta_g_paratec
        data = source.get('delta_g_total_paratec', None)
        if data is not None:
            delta_g_total_paratec = []
            for atom_data in data:
                values = data
                delta_g_total_paratec.append(values)
            out['delta_g_paratec'] = [dict(value=item) for item in delta_g_total_paratec]

        # delta_g
        data = source.get('delta_g_total', None)
        if data is not None:
            delta_g_total = []
            for atom_data in data:
                values = data
                delta_g_total.append(values)
            out['delta_g'] = [dict(value=item) for item in delta_g_total]


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
