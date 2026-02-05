import os
from pathlib import Path as PathLib
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

    def get_pseudopotentials_xml(
        self, arrays: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        """
        Idiomatic transformer: Extract limited pseudopotential metadata from
        vasprun.xml.

        Note: vasprun.xml contains very limited pseudopotential information
        compared to OUTCAR. Only name (TITEL) and valence electrons (ZVAL) are
        available. Other fields like type, XC functional, and cutoffs must be
        supplemented from OUTCAR via multi-pass parsing.

        Args:
            arrays: List of arrays from atominfo, need to filter for
                @name='atomtypes'

        Returns:
            list[dict]: List of pseudopotential dicts with limited metadata
                (name and n_valence_electrons only)
        """
        if not arrays:
            return []

        # Find the atomtypes array
        atomtypes_array = None
        for arr in arrays:
            if (
                isinstance(arr, dict)
                and arr.get(f'{self.attribute_prefix}name') == 'atomtypes'
            ):
                atomtypes_array = arr
                break

        if not atomtypes_array:
            return []

        pseudopotentials = []

        # Extract atomtypes data - vasprun.xml stores as array of rc/c elements
        # Structure: array[@name='atomtypes'] -> set -> rc with c elements for
        # atomspertype, element, pseudopotential, valence
        for atomtype_set in [atomtypes_array]:
            if not isinstance(atomtype_set, dict):
                continue

            # Get the 'set' element which contains 'rc' elements for each atom type
            rc_elements = atomtype_set.get('set', {}).get('rc', [])
            if not isinstance(rc_elements, list):
                rc_elements = [rc_elements] if rc_elements else []

            for rc in rc_elements:
                if not isinstance(rc, dict):
                    continue

                # Each rc has 'c' elements: c[0]=atomspertype, c[1]=element,
                # c[2]=mass, c[3]=valence, c[4]=pseudopotential
                c_elements = rc.get('c', [])
                min_elements = 5  # Need indices 0-4 for all pseudopotential data
                if not isinstance(c_elements, list) or len(c_elements) < min_elements:
                    continue

                # Extract name and valence electrons from c elements
                idx_pp_name = 4  # Index for pseudopotential name
                idx_valence = 3  # Index for valence electrons
                pp_name = (
                    c_elements[idx_pp_name] if len(c_elements) > idx_pp_name else None
                )
                valence_str = (
                    c_elements[idx_valence] if len(c_elements) > idx_valence else None
                )

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

        # First pass: XML_KEY for basic structure
        data_parser.annotation_key = vasp.XML_KEY
        xml_parser.convert(data_parser)

        # Second pass: XML2_KEY for additional XML data
        data_parser.annotation_key = vasp.XML2_KEY
        xml_parser.convert(data_parser)

        # Third pass: OUTCAR_KEY to extend with OUTCAR data if available
        # This allows OUTCAR to supplement vasprun.xml pseudopotentials with
        # detailed metadata
        outcar_path = self._find_outcar()
        if outcar_path and os.path.exists(outcar_path):
            LOGGER.info(
                f'Found OUTCAR at {outcar_path}, extending vasprun.xml data '
                'with detailed pseudopotential metadata'
            )
            from nomad_simulation_parsers.parsers.vasp.outcar_parser import (
                OutcarParser,
                OutcarTextParser,
            )

            outcar_parser = OutcarParser()
            outcar_parser.text_parser = OutcarTextParser()
            outcar_parser.filepath = outcar_path

            data_parser.annotation_key = vasp.OUTCAR_KEY
            # Merge by index position: OUTCAR PP[0] extends XML PP[0], etc.
            # This preserves XML structure while adding OUTCAR's detailed metadata
            outcar_parser.convert(data_parser, update_mode='merge')

            outcar_parser.close()

        self.archive.data = data_parser.data_object

        # close file objects
        data_parser.close()
        xml_parser.close()

    def _find_outcar(self) -> str | None:
        """Find OUTCAR file in the same directory as vasprun.xml."""
        mainfile_dir = PathLib(self.mainfile).parent

        # Check for any file starting with 'outcar' (case-insensitive)
        # Catches: OUTCAR, outcar, OUTCAR.gz, outcar.bz2, etc.
        for file in mainfile_dir.iterdir():
            if file.name.lower().startswith('outcar'):
                return str(file)

        return None
