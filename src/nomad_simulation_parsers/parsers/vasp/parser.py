import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel import EntryArchive
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.parsing import MatchingParser

from nomad_simulation_parsers.parsers.phonopy.utils import find_phonopy_yaml

from .outcar_parser import OutcarArchiveWriter
from .xml_parser import XMLArchiveWriter


class VASPParser(MatchingParser):
    level = 2
    creates_children = True

    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ):
        matched = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
        if matched:
            yaml_file = find_phonopy_yaml(filename, 'vasp')
            if yaml_file:
                yaml_key = os.path.relpath(yaml_file, os.path.dirname(filename))
                return [f'phonon_workflow:{yaml_key}']

        return matched

    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = {},
    ) -> None:
        if 'outcar' in mainfile.lower():
            archive_writer = OutcarArchiveWriter()
        else:
            archive_writer = XMLArchiveWriter()
        archive_writer.write(mainfile, archive, logger, child_archives)
