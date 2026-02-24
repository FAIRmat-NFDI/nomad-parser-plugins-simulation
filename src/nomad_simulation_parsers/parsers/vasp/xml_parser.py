import os
from pathlib import Path as PathLib
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.parsing.file_parser import ArchiveWriter, Quantity, TextParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.vasp.outcar_parser import potcar_quantities
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

    def to_float(self, value: str) -> float | None:
        """Convert string value to float (field-level transformer)."""
        return float(value) if value is not None else None

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


class XMLArchiveWriter(ArchiveWriter):
    def write_to_archive(self) -> None:
        data_parser = VASPMetainfoParser()
        data_parser.data_object = Simulation()

        xml_parser = VasprunParser(filepath=self.mainfile)

        # First pass: XML_KEY for basic structure
        data_parser.annotation_key = vasp.XML_KEY
        xml_parser.convert(data_parser)

        # Second pass: KPOINTS_XML for KSpace numerical settings
        data_parser.annotation_key = vasp.KPOINTS_XML
        xml_parser.convert(data_parser)

        # Third pass: PP_XML for Pseudopotentials from XML
        data_parser.annotation_key = vasp.PP_XML
        xml_parser.convert(data_parser)

        # Fourth pass: XML2_KEY for additional XML data
        data_parser.annotation_key = vasp.XML2_KEY
        xml_parser.convert(data_parser)

        # Fifth pass: PP_OUT to extend with OUTCAR pseudopotential metadata
        # This allows OUTCAR to supplement vasprun.xml pseudopotentials with
        # detailed metadata (SHA256, LPAW, LULTRA, etc.)
        outcar_path = self._find_outcar()
        if outcar_path and os.path.exists(outcar_path):
            LOGGER.info(
                f'Found OUTCAR at {outcar_path}, extending vasprun.xml '
                'pseudopotentials with detailed metadata'
            )
            potcar_pattern = (
                r'POTCAR:([\s\S]+?VRHFIN[\s\S]+?)'
                r'(?=\s*POTCAR:|\s*local pseudopotential:|\Z)'
            )
            outcar_supplement_parser = TextParser(
                quantities=[
                    Quantity(
                        'pseudopotentials',
                        potcar_pattern,
                        repeats=True,
                        sub_parser=TextParser(quantities=potcar_quantities),
                    )
                ]
            )

            outcar_parser = MappingTextParser(filepath=outcar_path)
            outcar_parser.text_parser = outcar_supplement_parser

            data_parser.annotation_key = vasp.PP_OUT
            # Merge by index position: OUTCAR PP[0] extends XML PP[0], etc.
            # This preserves XML structure while adding OUTCAR's detailed metadata
            outcar_parser.convert(data_parser, update_mode='merge')

            # Clean up duplicate pseudopotentials created by type mismatch
            # When PP_OUT merges at index 0 but finds KSpace, creates new PP
            model_method = data_parser.data_object.model_method[0]
            seen_pp_names = set()
            deduplicated_ns = []
            for ns in model_method.numerical_settings:
                if ns.m_def.name == 'Pseudopotential':
                    pp_name = getattr(ns, 'name', None)
                    if pp_name and pp_name in seen_pp_names:
                        LOGGER.debug(f'Removed duplicate Pseudopotential: {pp_name}')
                        continue
                    if pp_name:
                        seen_pp_names.add(pp_name)
                deduplicated_ns.append(ns)

            # Replace with deduplicated list
            model_method.numerical_settings = deduplicated_ns

            outcar_parser.close()

        self.archive.data = data_parser.data_object

        # close file objects
        data_parser.close()
        xml_parser.close()

    def _find_outcar(self) -> str | None:
        """Find OUTCAR file in the same directory as vasprun.xml.

        Matches any file starting with 'outcar' (case-insensitive):
        OUTCAR, outcar, OUTCAR.gz, outcar.bz2, etc.
        """
        mainfile_dir = PathLib(self.mainfile).parent
        return next(
            (
                str(f)
                for f in mainfile_dir.iterdir()
                if f.name.lower().startswith('outcar')
            ),
            None,
        )
