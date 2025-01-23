from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )

from nomad.parsing.file_parser import ArchiveWriter
from nomad_simulations.schema_packages.workflow import SinglePoint


class WorkflowWriter(ArchiveWriter):
    """
    Base class for parsing the workflow section.

    """

    archives: dict[str, 'EntryArchive'] = {}

    def write_workflow(self):
        """
        Abstract method to write workflow section.
        """
        if self.archive.workflow2 is not None:
            return

        outputs = self.archive.outputs
        if len(outputs) == 1:
            self.archive.workflow2 = SinglePoint()

    def write_to_archive(self):
        super().write_to_archive()
        self.write_workflow()


class DFTGWWorkflowWriter(WorkflowWriter):
    """
    Base class for DFT + GW workflows.
    """

    def write_workflow(self):
        """
        Connect the DFT and GW archives in the workflow section.
        """
        # TODO fill workflow2 once gw workflow schema is defined in nomad_simulations
        dft_archive = self.archive.get('dft')  # noqa
        gw_archive = self.archive.get('gw')  # noqa
