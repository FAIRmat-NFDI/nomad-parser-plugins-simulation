from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from importlib import reload

from nomad.parsing import MatchingParser

from nomad_simulation_parsers.schema_packages import vasp
from nomad_simulation_parsers.schema_packages.utils import remove_mapping_annotations

from .outcar_parser import OutcarArchiveWriter
from .xml_parser import XMLArchiveWriter


class VASPParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = {},
    ) -> None:
        # reload schema to load vasp annotations
        reload(vasp)

        if 'outcar' in mainfile.lower():
            archive_writer = OutcarArchiveWriter()
        else:
            archive_writer = XMLArchiveWriter()
        archive_writer.write(mainfile, archive, logger, child_archives)

        # remove annotations
        # TODO cache? put in close context
        remove_mapping_annotations(vasp.general.Simulation.m_def)
