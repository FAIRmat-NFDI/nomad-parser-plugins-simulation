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

        # electric field gradient
        data = source.get('efg', [])
        electric_field_gradients = []
        for i, atom_data in enumerate(data):
            values = np.reshape(atom_data[2:], (3, 3))
            sec_efg = self.e_field_gradient_class(
                type="total", entity_ref=particle_state[i]
            )
            sec_efg.value = values
            electric_field_gradients.append(sec_efg)
        out['electric_field_gradients'] = [dict(value=e) for e in electric_field_gradients]

        txt = ''
        for key, value in out.items():
            txt += f'{key}          {value}\n'

        with open('ciao.txt', 'w') as f:
            f.write(txt)
        # debug(out)
        
        return [out]



class GIPAWMainfileXMLParser(MainfileXMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER
    
    def get_nmr_xml(self, pippo: dict[str, Any], **kwargs) -> Any:
        from devtools import debug
        debug("sto facendo qualcosa")
        debug(pippo)

        """
        # shielding_tensors
        if 'ms_list' not in self._results:
            st = self.fileparser.results._data['gpw:gipaw[0]']['output[0]']['shielding_tensors[0]']
            ms_list = []
            for key, value in st.items():
                if not isinstance(value, dict):
                    continue

                for atom in value['_data']:
                    atom_list = []
                    atom_list.append(atom['name'])
                    atom_list.append(int(atom['index']))
                    atom_list = atom_list + self.extract_floats_from_string(atom['atom'])
                    ms_list.append(atom_list)
            
            self._results['ms_list'] = ms_list"""
        
        # magnetic shieldings
        source = self.data['{http://www.quantum-espresso.org/ns/gpw/qes_gipaw_1.0}gipaw']['output']['shielding_tensors']['atom']
        debug(source)
        data = []
        for atom_source in source:
            atom_list = []
            atom_list.append(atom_source['@name'])
            atom_list.append(int(atom_source['@index']))
            atom_list = atom_list + atom_source['__value']
            data.append(atom_list)

        debug(data)
        magnetic_shieldings = []
        for atom_data in data:
            values = np.reshape(atom_data[2:], (3, 3))
            FACTOR = 1e-6
            magnetic_shieldings.append(values * FACTOR * ureg("dimensionless"))
        out = dict(magnetic_shieldings=[dict(value=m) for m in magnetic_shieldings ])
        debug(out)

        return [out]



class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    _text_parser = GIPAWMainfileTextParser(text_parser=GIPAWFileParser())
    _xml_parser = GIPAWMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
