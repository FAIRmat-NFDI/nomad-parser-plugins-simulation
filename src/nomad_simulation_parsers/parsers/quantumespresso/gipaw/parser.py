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

    def get_gipaw_text(self, source: dict[str, Any]) -> list[dict[str, Any]]:
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
        efg = source.get('efg', None)
        if efg is not None:
            electric_field_gradients = []
            for atom_data in efg:
                values = np.reshape(atom_data[2:], (3, 3))
                electric_field_gradients.append(values)
            out['electric_field_gradients'] = [dict(value=item) for item in electric_field_gradients]

        # hyperfine_dipolar
        hd = source.get('hyperfine_dipolar', None)
        if hd is not None:
            hyperfine_dipolar = []
            for atom_data in hd:
                values = np.reshape(atom_data[2:], (3, 3))
                hyperfine_dipolar.append(values)
            out['hyperfine_dipolar'] = [dict(value=item) for item in hyperfine_dipolar]

        # hyperfine_fermi_contact
        hfc = source.get('hyperfine_fermi_contact', None)
        if hfc is not None:
            hyperfine_fermi_contact = []
            for atom_data in hfc:
                values = atom_data[-1]
                hyperfine_fermi_contact.append(values)
            out['hyperfine_fermi_contact'] = [dict(value=item) for item in hyperfine_fermi_contact]

        # delta_g_paratec
        delta_g_paratec = source.get('delta_g_total_paratec', None)
        if delta_g_paratec is not None:
            out['delta_g_paratec'] = dict(value=delta_g_paratec)

        # delta_g
        delta_g = source.get('delta_g_total', None)
        if delta_g is not None:
            out['delta_g'] = dict(value=delta_g)

        return [out]



class GIPAWMainfileXMLParser(MainfileXMLParser):
    _job: str | None = None

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    @property
    def job(self) -> str | None:
        if self._job is None:
            try:
                self._job = self.data["input"]["job"]
            except Exception as exc:
                self.logger.warning("Unable to get job from data: %s", exc)
        return self._job


    def get_magnetic_shieldings(self, atom: dict[str, Any]) -> Any:
        if self.job != "nmr":
            return
        value = np.reshape(atom.get("__value"), (3, 3))
        FACTOR = 1e-6
        return value * FACTOR * ureg("dimensionless")

    def get_magnetic_susceptibilities(self, source: dict[str, Any], **kwargs) -> Any:
        if self.job != "nmr":
            return
        name = kwargs.get("name")
        if name != "value":
            value = source.get(name, None)
            return np.reshape(value.get("__value", None), (3, 3))

        value_vgv = source.get("susceptibility_low", None)
        vgv = np.reshape(value_vgv.get("__value", None), (3, 3))
        value_pgv = source.get("susceptibility_high", None)
        pgv = np.reshape(value_pgv.get("__value", None), (3, 3))
        sus = (vgv + pgv) / 2
        return sus

    def get_efg(self, atom: dict[str, Any]) -> Any:
        if self.job != "efg":
            return
        return np.reshape(atom.get("__value"), (3, 3))

    def get_hyperfine_dipolar(self, atom: dict[str, Any]) -> Any:
        if self.job != "hyperfine":
            return
        return np.reshape(atom.get("__value"), (3, 3))

    def get_hyperfine_fermi_contact(self, atom: dict[str, Any]) -> Any:
        if self.job != "hyperfine":
            return
        return atom.get("__value")

    def get_delta_g(self, source: dict[str, Any]) -> Any:
        if self.job != "g-tensor":
            return
        return np.reshape(source.get("__value"), (3, 3))




class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    _text_parser = GIPAWMainfileTextParser(text_parser=GIPAWFileParser())
    _xml_parser = GIPAWMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
