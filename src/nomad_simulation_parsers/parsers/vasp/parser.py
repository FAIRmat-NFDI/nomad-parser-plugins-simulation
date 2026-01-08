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

from nomad_simulation_parsers.schema_packages import (
    vasp,  # noqa: F401 - needed to register mapping annotations
)

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
        # TEMPORARY: reload to work around annotation contamination
        from nomad_simulations.schema_packages import general, model_method, numerical_settings

        print(f"\nVASP Before reload:")
        print(f"  general.Simulation.m_def: {general.Simulation.m_def.m_annotations.get('mapping', {}).keys()}")
        print(f"  model_method.DFT.xc: {model_method.DFT.xc.m_annotations.get('mapping', {}).keys()}")

        # Check Pseudopotential class BEFORE reload
        pseudopot_class_before = vasp.Pseudopotential
        pseudopot_mdef_before = vasp.Pseudopotential.m_def
        pseudopot_annotations_before = vasp.Pseudopotential.m_def.m_annotations.get('mapping', {})
        print(f"  vasp.Pseudopotential.m_def (BEFORE): {pseudopot_annotations_before.keys()}")
        print(f"  vasp.Pseudopotential class id (BEFORE): {id(pseudopot_class_before)}")
        print(f"  vasp.Pseudopotential.m_def id (BEFORE): {id(pseudopot_mdef_before)}")

        # Check what NOMAD's plugin system references
        print(f"  numerical_settings.Pseudopotential class id: {id(numerical_settings.Pseudopotential)}")
        print(f"  numerical_settings.Pseudopotential.m_def id: {id(numerical_settings.Pseudopotential.m_def)}")
        print(f"  numerical_settings.Pseudopotential annotations: {numerical_settings.Pseudopotential.m_def.m_annotations.get('mapping', {}).keys()}")

        reload(vasp)

        print(f"\nVASP After reload:")
        print(f"  general.Simulation.m_def: {general.Simulation.m_def.m_annotations.get('mapping', {}).keys()}")
        print(f"  model_method.DFT.xc: {model_method.DFT.xc.m_annotations.get('mapping', {}).keys()}")

        # Check Pseudopotential class AFTER reload
        pseudopot_class_after = vasp.Pseudopotential
        pseudopot_mdef_after = vasp.Pseudopotential.m_def
        pseudopot_annotations_after = vasp.Pseudopotential.m_def.m_annotations.get('mapping', {})
        print(f"  vasp.Pseudopotential.m_def (AFTER): {pseudopot_annotations_after.keys()}")
        print(f"  vasp.Pseudopotential class id (AFTER): {id(pseudopot_class_after)}")
        print(f"  vasp.Pseudopotential.m_def id (AFTER): {id(pseudopot_mdef_after)}")

        # Check if NOMAD still references the OLD class
        print(f"  numerical_settings.Pseudopotential class id (still): {id(numerical_settings.Pseudopotential)}")
        print(f"  numerical_settings.Pseudopotential.m_def id (still): {id(numerical_settings.Pseudopotential.m_def)}")
        print(f"  numerical_settings.Pseudopotential annotations (still): {numerical_settings.Pseudopotential.m_def.m_annotations.get('mapping', {}).keys()}")

        print(f"\n  Class changed? {pseudopot_class_before is not pseudopot_class_after}")
        print(f"  m_def changed? {pseudopot_mdef_before is not pseudopot_mdef_after}")
        print(f"  NOMAD references old class? {id(numerical_settings.Pseudopotential) == id(pseudopot_class_before)}")
        print(f"  NOMAD references new class? {id(numerical_settings.Pseudopotential) == id(pseudopot_class_after)}\n")

        if 'outcar' in mainfile.lower():
            archive_writer = OutcarArchiveWriter()
        else:
            archive_writer = XMLArchiveWriter()
        archive_writer.write(mainfile, archive, logger, child_archives)

        # NOTE: Removed remove_mapping_annotations() - Since we're no longer reloading,
        # the annotations are loaded once at module import and should persist across
        # all parses. Removing them would break subsequent parses.
