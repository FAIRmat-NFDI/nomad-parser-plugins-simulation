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

from nomad_simulation_parsers.parsers.utils.general import remove_mapping_annotations
from nomad_simulation_parsers.schema_packages import vasp

from .outcar_parser import OutcarArchiveWriter
from .xml_parser import XMLArchiveWriter


def ref_reciprocal_lattice(archive: 'EntryArchive', logger: 'BoundLogger') -> None:
    try:
        d = archive.data
        recip = d.model_method[0].numerical_settings[0].reciprocal_lattice_vectors
        if recip is not None and len(recip) > 0:
            d.outputs[0].electronic_eigenvalues[0].reciprocal_cell = recip
    except (AttributeError, IndexError) as e:
        logger.warning(f'Failed to set reciprocal lattice vectors: {e}')


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

        # ref_reciprocal_lattice(archive, logger)

        # remove annotations
        # TODO cache? put in close context
        remove_mapping_annotations(vasp.general.Simulation.m_def)
