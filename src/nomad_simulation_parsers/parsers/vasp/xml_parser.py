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

    def get_pseudopotentials(self, atomtypes: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """
        Extract basic pseudopotential data from vasprun.xml atomtypes array.

        The atomtypes array contains:
        - atomspertype (count)
        - element
        - mass
        - valence (n_valence_electrons)
        - pseudopotential (name/TITEL)

        Note: vasprun.xml does NOT contain detailed POTCAR metadata (LPAW, LULTRA,
        LEXCH, ENMAX, ENMIN, RCORE, VRHFIN, SHA256). These are only available in
        OUTCAR or by parsing POTCAR directly.
        """
        if not atomtypes:
            LOGGER.debug('get_pseudopotentials: No atomtypes provided')
            return []

        LOGGER.debug(f'get_pseudopotentials: Processing {len(atomtypes)} atomtypes')

        pseudopotentials = []
        for atomtype in atomtypes:
            # Extract fields from the rc array structure
            # Format: [atomspertype, element, mass, valence, pseudopotential_name]
            rc = atomtype.get('rc', [])
            if not rc or len(rc) < 5:
                LOGGER.debug(f'get_pseudopotentials: Skipping atomtype with rc={rc}')
                continue

            pp_data = {
                'name': rc[4].strip(),  # pseudopotential name (TITEL)
                'n_valence_electrons': float(rc[3]),  # valence
            }
            pseudopotentials.append(pp_data)
            LOGGER.debug(f'get_pseudopotentials: Added {pp_data["name"]}')

        LOGGER.debug(f'get_pseudopotentials: Returning {len(pseudopotentials)} pseudopotentials')
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

        # TODO: Add POTCAR parser for complete pseudopotential metadata
        # Currently vasprun.xml provides only basic info (name, valence) while OUTCAR
        # contains full POTCAR headers. Direct POTCAR parsing would enable complete
        # pseudopotential support (LPAW, LULTRA, LEXCH, ENMAX/ENMIN, RCORE, VRHFIN,
        # SHA256) regardless of which mainfile is used. This would allow proper type
        # determination and XC functional resolution from vasprun.xml-only uploads.

        # close file objects
        data_parser.close()
        xml_parser.close()
