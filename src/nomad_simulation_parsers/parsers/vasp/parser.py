from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad_simulation_parsers.parsers.utils import remove_mapping_annotations
from .outcar_parser import VASPOutcarParser
from .xml_parser import VASPXMLParser


class VASPParser:
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:

        # import schema to load annotations
        from nomad_simulation_parsers.schema_packages import vasp

        if 'outcar' in mainfile.lower():
            parser = VASPOutcarParser()
        else:
            parser = VASPXMLParser()
        # TODO remove this for debug
        self.parser = parser
        parser.parse(mainfile, archive, logger, child_archives)

        # remove annotations
        # TODO cache? put in close context
        remove_mapping_annotations(vasp.general.Simulation.m_def)
