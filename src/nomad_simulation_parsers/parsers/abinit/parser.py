from collections.abc import Iterable
from importlib import reload
from typing import Any

from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import (
    DFTGWWorkflow,
    MolecularDynamics,
    SinglePoint,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
    GeometryOptimizationModel,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import abinit

from .file_parser import AbinitOutParser

LOGGER = get_logger(__name__)


# TODO temporary fix for structlog unable to propagate logger
class AbinitMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    text_parser = AbinitOutParser()

    def get_workflow_method(self) -> str:
        ionmov = self.data_object.get_input_var('ionmov', 1, [0])[0]
        return {
            1: 'viscous_damped_md',
            2: 'bfgs',
            3: 'bfgs',
            4: 'conjugate_gradient',
            5: 'steepest_descent',
            7: 'quenched_md',
            10: 'dic_bfgs',
            11: 'dic_bfgs',
            20: 'diis',
        }.get(ionmov)

    def get_input_var(self, **kwargs: dict[str, Any]):
        return self.data_object.get_input_var(
            kwargs.get('name'), kwargs.get('n_dataset', 1), kwargs.get('default')
        )


class AbinitArchiveWriter(ArchiveWriter):
    mainfile_parser = MainfileParser()
    metainfo_parser = MetainfoParser()
    code_name = 'ABINIT'
    annotation_key = 'out'

    def parse_workflow(self):
        ionmov = self.mainfile_parser.data_object.get_input_var('ionmov', 1, [0])[0]
        vis = self.mainfile_parser.data_object.get_input_var('vis', 1, [100.0])[0]
        if ionmov in [2, 3, 4, 5, 7, 10, 11, 20] or (ionmov == 1 and vis > 0.0):
            self.archive.workflow2 = GeometryOptimization(
                model=GeometryOptimizationModel()
            )
        elif ionmov in [6, 8, 9, 12, 13, 14, 23] or (ionmov == 1 and vis == 0.0):
            self.archive.workflow2 = MolecularDynamics()
        else:
            self.archive.workflow2 = SinglePoint()
        self.metainfo_parser.annotation_key = self.annotation_key
        self.metainfo_parser.data_object = self.archive.workflow2
        self.mainfile_parser.convert(self.metainfo_parser)

    def write_to_archive(self):
        reload(abinit)

        self.archive.data = Simulation(program=Program(name=self.code_name))
        self.metainfo_parser.annotation_key = self.annotation_key
        self.metainfo_parser.data_object = self.archive.data

        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.convert(self.metainfo_parser)

        self.parse_workflow()

        gw_archive = self.child_archives.get('GW')
        if gw_archive is not None:
            gw_archive.data = Simulation(program=Program(name=self.code_name))

            writer = AbinitArchiveWriter()
            writer.annotation_key = 'gw_out'
            writer.write(self.mainfile, gw_archive, self.logger)

            workflow_archive = self.child_archives['GW_workflow']
            workflow_archive.workflow2 = DFTGWWorkflow(
                tasks=[self.archive.workflow2, gw_archive.workflow2]
            )


class AbinitParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = AbinitArchiveWriter()

    def is_mainfile(
        self,
        filename: str,
        mime,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable:
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
        if is_mainfile:
            out_parser = AbinitOutParser()
            out_parser.findall = False
            out_parser.mainfile = filename
            ds_numbers = out_parser.dataset_numbers
            optdriver = out_parser.input_vars.get('optdriver', [])
            out_parser.findall = True
            n_gw = [4, 66]
            if n_gw[0] in ds_numbers and (1 and 2 and 3) not in ds_numbers:
                return True
            if len(optdriver) == n_gw[0] and (optdriver[-1] in n_gw):
                self.creates_children = True
                return ['GW', 'GW_workflow']
            return True

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
