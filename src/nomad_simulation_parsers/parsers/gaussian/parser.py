from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from importlib import reload

from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.schema_packages import gaussian

from .text_parser import GaussianOutReader

LOGGER = get_logger(__name__)


class OutParser(MappingTextParser):
    """
    Couples GaussianOutReader (regex) with a few convenience getters that
    the mapping rules will call.
    """

    def __init__(self):
        super().__init__(text_parser=GaussianOutReader())

    def get_program_data(self, src: dict[str, Any]) -> dict[str, Any]:
        return {
            'program_name': 'Gaussian',
            #'program_version': src.get('program_version'),
        }


class GaussianParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        reload(gaussian)

        reader = OutParser()
        reader.filepath = mainfile

        meta = MetainfoParser(data_object=Simulation())
        meta.annotation_key = 'out'
        # meta.max_nested_level = 1

        reader.convert(meta)
        archive.data = meta.data_object

        meta.close()
        reader.close()
