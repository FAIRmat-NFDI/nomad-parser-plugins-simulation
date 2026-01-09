from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.parsing import MatchingParser

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
        # NOTE: Removed reload(vasp) - it was breaking metainfo registration by creating
        # new class instances with different m_def objects. The annotations are loaded
        # when the module is first imported and don't need to be reloaded.

        if 'outcar' in mainfile.lower():
            archive_writer = OutcarArchiveWriter()
        else:
            archive_writer = XMLArchiveWriter()
        archive_writer.write(mainfile, archive, logger, child_archives)

        # NOTE: Removed remove_mapping_annotations() - Since we're no longer reloading,
        # the annotations are loaded once at module import and should persist across
        # all parses. Removing them would break subsequent parses.
