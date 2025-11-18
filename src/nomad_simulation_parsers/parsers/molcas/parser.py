from importlib import reload
from typing import TYPE_CHECKING, Any

from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.schema_packages import molcas

from .text_parser import MolcasOutReader

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import EntryArchive
    from structlog.stdlib import BoundLogger

LOGGER = get_logger(__name__)


class OutParser(MappingTextParser):
    """
    Couples MolcasOutReader with convenience getters used in the mapping rules.
    """

    def __init__(self) -> None:
        super().__init__(text_parser=MolcasOutReader())

    def get_program_data(self, src: dict[str, Any]) -> dict[str, Any]:
        return {
            'program_name': 'OpenMolcas',
        }


class MolcasParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        # Hot reload
        reload(molcas)

        reader = OutParser()
        reader.filepath = mainfile

        meta = MetainfoParser(data_object=Simulation())
        meta.annotation_key = 'out'

        reader.convert(meta)
        archive.data = meta.data_object

        meta.close()
        reader.close()
