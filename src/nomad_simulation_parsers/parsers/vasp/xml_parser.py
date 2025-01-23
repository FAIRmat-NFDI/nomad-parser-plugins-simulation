from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.utils.general import remove_mapping_annotations


class VasprunParser(XMLParser):
    def mix_alpha(self, mix: float, cond: bool) -> float:
        return mix if cond else 0

    def get_eigenvalues(self, array: list) -> dict[str, Any]:
        if array is None:
            return {}
        transposed = np.transpose(array)
        return dict(eigenvalues=transposed[0], occupations=transposed[1])

    def get_energy_contributions(
        self, source: dict[str, Any], **kwargs
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


class XMLArchiveWriter(ArchiveWriter):
    def write_to_archive(self) -> None:
        # import schema to load annotations
        from nomad_simulation_parsers.schema_packages import vasp

        data_parser = MetainfoParser()
        data_parser.data_object = Simulation()

        xml_parser = VasprunParser(filepath=self.mainfile)

        data_parser.annotation_key = 'xml'
        xml_parser.convert(data_parser)

        data_parser.annotation_key = 'xml2'
        xml_parser.convert(data_parser)

        self.archive.data = data_parser.data_object

        # close file objects
        data_parser.close()
        xml_parser.close()

        # remove annotations
        # TODO cache? put in close context
        remove_mapping_annotations(vasp.general.Simulation.m_def)
