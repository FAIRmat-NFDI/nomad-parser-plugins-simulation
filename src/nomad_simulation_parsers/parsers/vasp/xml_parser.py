from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.schema_packages import vasp

LOGGER = get_logger(__name__)


# TODO temporary fix for structlog unable to propagate logger
class VASPMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class VasprunParser(XMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def mix_alpha(self, mix: float, cond: bool) -> float:
        return mix if cond else 0

    def get_eigenvalues(self, array: list) -> dict[str, Any]:
        if array is None:
            return {}
        transposed = np.transpose(array)
        return dict(eigenvalues=transposed[0], occupations=transposed[1])

    def get_energy_contributions(
        self, source: list[dict[str, Any]], **kwargs
    ) -> list[dict[str, Any]]:
        return [
            c
            for c in source
            if c.get(f'{self.attribute_prefix}name') not in kwargs.get('exclude', [])
        ]

    def get_data(self, source: dict[str, Any], **kwargs) -> Any:
        if source.get(self.value_key):
            return source[self.value_key]
        path = kwargs.get('path')
        if path is None:
            return

        parser = Path(path=path)
        return parser.get_data(source)

    def get_forces(self, source: dict[str, Any]) -> dict[str, Any]:
        value = self.get_data(source, path='.varray.v')
        if value is None:
            return {}
        return dict(forces=value, npoints=len(value), rank=[3])

    def reshape_array(self, source: np.ndarray, shape_rest: tuple = (3,)) -> np.ndarray:
        if source is None:
            return
        return np.reshape(
            source, (np.size(source) // int(np.prod(shape_rest)), *shape_rest)
        )

    def get_pseudopotentials_xml(self, atomtypes_array: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """
        Idiomatic transformer: Extract limited pseudopotential metadata from vasprun.xml.

        Note: vasprun.xml contains very limited pseudopotential information compared to OUTCAR.
        Only name (TITEL) and valence electrons (ZVAL) are available. Other fields like type,
        XC functional, and cutoffs must be supplemented from OUTCAR via multi-pass parsing.

        Args:
            atomtypes_array: The atominfo.array with @name='atomtypes' from vasprun.xml

        Returns:
            list[dict]: List of pseudopotential dicts with limited metadata (name and n_valence_electrons only)
        """
        if not atomtypes_array:
            return []

        pseudopotentials = []

        # Extract atomtypes data - vasprun.xml stores as array of rc/c elements
        # Structure: array[@name='atomtypes'] -> set -> rc with c elements for atomspertype, element, pseudopotential, valence
        for atomtype_set in atomtypes_array:
            if not isinstance(atomtype_set, dict):
                continue

            # Get the 'set' element which contains 'rc' elements for each atom type
            rc_elements = atomtype_set.get('set', {}).get('rc', [])
            if not isinstance(rc_elements, list):
                rc_elements = [rc_elements] if rc_elements else []

            for rc in rc_elements:
                if not isinstance(rc, dict):
                    continue

                # Each rc has 'c' elements: c[0]=atomspertype, c[1]=element, c[2]=pseudopotential name, c[3]=valence
                c_elements = rc.get('c', [])
                if not isinstance(c_elements, list) or len(c_elements) < 4:
                    continue

                # Extract name and valence electrons from c elements
                pp_name = c_elements[2] if len(c_elements) > 2 else None
                valence_str = c_elements[3] if len(c_elements) > 3 else None

                # Parse valence electrons
                n_valence = None
                if valence_str is not None:
                    try:
                        n_valence = float(valence_str)
                    except (ValueError, TypeError):
                        pass

                # Create pseudopotential dict with limited metadata
                pp_data = {
                    'name': pp_name,
                    'n_valence_electrons': n_valence,
                }

                pseudopotentials.append(pp_data)

        return pseudopotentials


class XMLArchiveWriter(ArchiveWriter):
    def write_to_archive(self) -> None:
        data_parser = VASPMetainfoParser()
        data_parser.data_object = Simulation()

        xml_parser = VasprunParser(filepath=self.mainfile)

        data_parser.annotation_key = vasp.XML_KEY
        xml_parser.convert(data_parser)

        data_parser.annotation_key = vasp.XML2_KEY
        xml_parser.convert(data_parser)

        self.archive.data = data_parser.data_object

        # close file objects
        data_parser.close()
        xml_parser.close()
