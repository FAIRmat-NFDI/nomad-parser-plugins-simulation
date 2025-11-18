from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import EntryArchive
    from structlog.stdlib import BoundLogger

from importlib import reload

from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.schema_packages import pyscf

from .text_parser import PySCFOutReader

LOGGER = get_logger(__name__)


class OutParser(MappingTextParser):
    """
    Couples PySCFOutReader (regex) with a few convenience getters that
    the mapping rules will call.
    """

    def __init__(self) -> None:
        super().__init__(text_parser=PySCFOutReader())

    def get_program_data(self, src: dict[str, Any]) -> dict[str, Any]:
        return {
            'program_name': 'PySCF',
            # 'program_version': src.get('program_version'),
        }


class PySCFParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        # reload mapping rules on every parse, just like in gaussian.py
        reload(pyscf)

        reader = OutParser()
        reader.filepath = mainfile

        meta = MetainfoParser(data_object=Simulation())
        meta.annotation_key = 'out'  # same annotation key convention as Gaussian

        reader.convert(meta)
        archive.data = meta.data_object

        meta.close()
        reader.close()
