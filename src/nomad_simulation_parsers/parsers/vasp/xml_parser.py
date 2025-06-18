from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

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

    def get_eigenvalues(self, package: dict) -> dict[str, Any]:
        """
        Extracts eigenvalues and occupations from the VASP XML <eigenvalues.array> branch.
        """
        k_dict = next(
            filter(
                lambda x: x.get('__value', '') == 'kpoint',
                package.get('dimension', []),
            ),
            {},
        )
        k_level = int(k_dict.get('@dim', '0'))

        layer = package.get('set', {})
        for level in range(0, 3):
            if k_level == level:
                break
            else:
                layer = layer.get('set', {})

        # TODO: handle more lower layers
        data = np.transpose([lyr.get('r', []) for lyr in layer])
        return dict(eigenvalues=data[0], occupations=data[1])

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
        return dict(forces=value, npoints=len(value))  # ! remove npoints

    def reshape_array(self, source: np.ndarray, shape_rest: tuple = (3,)) -> np.ndarray:
        if source is None:
            return
        return np.reshape(
            source, (np.size(source) // int(np.prod(shape_rest)), *shape_rest)
        )

    def get_dos(self, source: list[list[float]] | None) -> dict[str, Any]:
        if source is None:
            return {}
        source = np.transpose(source)
        return dict(energies=source[0], value=source[1])


class XMLArchiveWriter(ArchiveWriter):
    def write_to_archive(self) -> None:
        data_parser = VASPMetainfoParser()
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
