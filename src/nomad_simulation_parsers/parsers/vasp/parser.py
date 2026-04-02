import os
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
    def _has_electronic_outputs(self, archive: 'EntryArchive') -> bool:
        outputs = getattr(getattr(archive, 'data', None), 'outputs', None) or []
        for output in outputs:
            if getattr(output, 'electronic_band_structures', None):
                return True
            if getattr(output, 'electronic_dos', None):
                return True
            if getattr(output, 'electronic_band_gaps', None):
                return True
        return False

    def _backfill_from_outcar(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'],
    ) -> None:
        if self._has_electronic_outputs(archive):
            return

        outcar_path = os.path.join(os.path.dirname(mainfile), 'OUTCAR')
        if not os.path.isfile(outcar_path):
            return

        from nomad.datamodel import EntryArchive

        outcar_archive = EntryArchive()
        OutcarArchiveWriter().write(outcar_path, outcar_archive, logger, child_archives)

        outcar_outputs = getattr(getattr(outcar_archive, 'data', None), 'outputs', None) or []
        if not outcar_outputs:
            return

        outcar_output = outcar_outputs[0]
        target_outputs = getattr(getattr(archive, 'data', None), 'outputs', None) or []
        if not target_outputs:
            archive.data.outputs = [outcar_output]
            return

        target_output = target_outputs[0]
        for quantity_name in (
            'electronic_band_structures',
            'electronic_band_gaps',
            'electronic_dos',
        ):
            if getattr(target_output, quantity_name, None):
                continue
            value = getattr(outcar_output, quantity_name, None)
            if value:
                setattr(target_output, quantity_name, value)

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

        # Some VASP XML datasets do not expose electronic outputs through current XML
        # extraction while equivalent OUTCAR payload is available in the same folder.
        if 'outcar' not in mainfile.lower():
            self._backfill_from_outcar(mainfile, archive, logger, child_archives)
