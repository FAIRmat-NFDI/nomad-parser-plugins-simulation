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
    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool:
        """
        Override mainfile detection to recognize standalone OUTCAR files.

        This override is necessary because:
        1. OUTCAR files don't match the XML pattern in the entry point regex
        2. We need content validation for standalone OUTCAR (not just filename)
        3. Works in conjunction with mainfile_alternative=True for auxiliary files

        Two scenarios:
        - Standalone OUTCAR: This method validates by content (vasp signature)
        - Auxiliary OUTCAR: Discovered via _find_outcar() in xml_parser.py

        Args:
            filename: Name of the file being checked
            mime: MIME type of the file
            buffer: Raw file contents as bytes
            decoded_buffer: Decoded file contents as string
            compression: Compression format if any

        Returns:
            True if this is a valid VASP mainfile (OUTCAR or vasprun.xml)
        """
        # Recognize standalone OUTCAR files by checking filename and content
        if 'OUTCAR' in filename.upper():
            # Validate it's a VASP OUTCAR by checking for VASP signatures
            # OUTCAR files contain lines like "vasp.X.X.X" and "executed on"
            buffer_lower = buffer.lower()
            return b'vasp' in buffer_lower and b'executed on' in buffer_lower

        # For vasprun.xml and other files, use default matching from entry point
        return super().is_mainfile(filename, mime, buffer, decoded_buffer, compression)

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
