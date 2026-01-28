import os
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.parsers.vasp.outcar_parser import (
    OutcarArchiveWriter,
    OutcarParser,
    OutcarTextParser,
)
from nomad_simulation_parsers.schema_packages import vasp

LOGGER = get_logger(__name__)

# Number of expected elements in atomtype rc data:
# [atomspertype, element, mass, valence, pseudopotential_name]
ATOMTYPE_RC_EXPECTED_LENGTH = 5


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

    def normalize_xc_label(self, raw_label: str) -> str:
        """
        Normalize VASP XC codes to standard functional names.

        Args:
            raw_label: VASP LEXCH code (e.g., 'PE' for PBE)

        Returns:
            Normalized functional name (e.g., 'GGA_X_PBE')
        """
        # Map VASP codes to libxc-style functional names
        vasp_xc_map = {
            'PE': 'GGA_X_PBE',  # PBE GGA
            'PS': 'GGA_X_PBE_SOL',  # PBEsol
            'CA': 'LDA_C_PZ',  # Perdew-Zunger LDA
            '91': 'GGA_X_PW91',  # Perdew-Wang 91
            'AM': 'GGA_X_AM05',  # AM05
            'RP': 'GGA_X_RPBE',  # RPBE
            'PW': 'GGA_X_PW91',  # PW91 (same as 91)
        }
        return vasp_xc_map.get(raw_label, raw_label)

    def get_pseudopotentials_xml(
        self, atomtypes: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        """
        Extract basic pseudopotential data from vasprun.xml atomtypes array.

        Args:
            atomtypes: List containing atomtypes array dict from vasprun.xml.
                      JMESPath filters to array with @name=='atomtypes'.

        Returns:
            List of dicts with keys: 'name', 'n_valence_electrons'

        Note: vasprun.xml does NOT contain detailed POTCAR metadata (LPAW, LULTRA,
        LEXCH, ENMAX, ENMIN, RCORE, VRHFIN, SHA256). These are only available in
        OUTCAR or by parsing POTCAR directly.
        """
        if not atomtypes or len(atomtypes) == 0:
            return []

        # Extract first array element (atomtypes is list with single dict)
        atomtype_array = atomtypes[0] if isinstance(atomtypes, list) else atomtypes
        rc_rows = atomtype_array.get('set', {}).get('rc', [])
        if not isinstance(rc_rows, list):
            rc_rows = [rc_rows]

        # Validate structure and extract pseudopotential data
        pseudopotentials = []
        for rc in rc_rows:
            c_elements = rc.get('c', [])
            if len(c_elements) >= ATOMTYPE_RC_EXPECTED_LENGTH:
                # Format: [atomspertype, element, mass, valence, pseudopotential_name]
                pseudopotentials.append(
                    {
                        'name': c_elements[4].strip(),
                        'n_valence_electrons': float(c_elements[3]),
                    }
                )

        return pseudopotentials


class XMLArchiveWriter(ArchiveWriter):
    def _supplement_pseudopotentials_from_outcar(
        self,
        archive_data: Simulation,
        outcar_parser_data: dict[str, Any],
    ) -> None:
        """
        Post-process pseudopotentials after OUTCAR supplementing.

        This method delegates to OutcarArchiveWriter._process_pseudopotentials
        to handle complex transformations:
        - Add PPCutoff from ENMAX/ENMIN
        - Determine type from LPAW/LULTRA flags
        - Map LEXCH codes to standard XC functionals
        - Link pseudopotentials to AtomsState

        Args:
            archive_data: The simulation archive being populated
            outcar_parser_data: Raw parsed OUTCAR data containing lpaw, lultra, lexch flags
        """
        # Delegate to OutcarArchiveWriter's post-processing logic
        # This ensures consistent behavior between standalone OUTCAR parsing
        # and XML+OUTCAR supplementing
        outcar_writer = OutcarArchiveWriter()
        outcar_writer._process_pseudopotentials(archive_data, outcar_parser_data)

    def write_to_archive(self) -> None:
        data_parser = VASPMetainfoParser()
        data_parser.data_object = Simulation()

        xml_parser = VasprunParser(filepath=self.mainfile)

        # DFT-specific key first to create DFT object
        data_parser.annotation_key = vasp.DFT_XML_KEY
        xml_parser.convert(data_parser)

        # Parse OUTCAR FIRST to create pseudopotentials in empty numerical_settings
        # This avoids positional merge bug where OUTCAR data would merge into XML's KSpace
        maindir = os.path.dirname(self.mainfile)
        outcar_files = search_files('OUTCAR', maindir)

        if outcar_files and data_parser.data_object.model_method:
            LOGGER.info(f'Parsing OUTCAR auxiliary file for pseudopotentials: {outcar_files[0]}')

            outcar_parser = OutcarParser()
            outcar_parser.text_parser = OutcarTextParser()
            outcar_parser.filepath = outcar_files[0]

            # Parse with default update_mode='merge' (PPs go into empty list)
            data_parser.annotation_key = vasp.OUTCAR_KEY
            outcar_parser.convert(data_parser)

            # Post-process pseudopotentials to add derived fields
            # (type, XC functional, cutoffs, atom linking)
            self._supplement_pseudopotentials_from_outcar(
                data_parser.data_object, outcar_parser.data
            )

            outcar_parser.close()

        # Parse XML_KEY AFTER OUTCAR
        # NOTE: XML_KEY collection annotations (KSpace) will not be added due to positional
        # merge bug. KSpace creation would need manual handling if required.
        data_parser.annotation_key = vasp.XML_KEY
        xml_parser.convert(data_parser)

        data_parser.annotation_key = vasp.XML2_KEY
        xml_parser.convert(data_parser)

        self.archive.data = data_parser.data_object

        # close file objects
        data_parser.close()
        xml_parser.close()
