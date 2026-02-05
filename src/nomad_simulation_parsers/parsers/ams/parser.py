import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import ams

from .file_parser import OutParser
from .file_parser import RKFParser as RKFTextParser

LOGGER = get_logger(__name__)


class MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_xc_functionals(self, source: dict[str, Any]) -> list[str]:
        xc_functionals = []
        for xc_type in ['LDA', 'GGA', 'MGGA']:
            functionals = source.get(xc_type, '').split()
            kind = ['XC'] if len(functionals) == 1 else ['X', 'C']
            for n, functional in enumerate(functionals):
                xc_functionals.append(
                    f'{xc_type}_{kind[n]}_{functional.rstrip("x").rstrip("c").upper()}'
                )
        return xc_functionals

    def get_contributions(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        return [
            dict(name=key, value=val) for key, val in source.items() if val is not None
        ]

    def get_eigenvalues(
        self, source: dict[str, Any] | list[np.ndarray]
    ) -> list[dict[str, np.ndarray]]:
        energies = source.get('energies', []) if isinstance(source, dict) else []
        occupations = (
            source.get('occupations', []) if isinstance(source, dict) else source[2]
        )
        nspin = max(len(energies), len(occupations))
        eigenvalues = [dict() for _ in range(nspin)]
        for n, energy in enumerate(energies):
            eigenvalues[n]['eigenvalues'] = energy
        for n, occupations in enumerate(occupations):
            eigenvalues[n]['occupations'] = occupations
        return [eig for eig in eigenvalues if eig]


class RKFParser(MainfileParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


# TODO temporary fix for structlog unable to propagate logger
class AMSMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class AMSArchiveWriter(ArchiveWriter):
    mainfile_parser = MainfileParser(text_parser=OutParser())
    metainfo_parser = AMSMetainfoParser()
    rkf_parser = RKFParser(text_parser=RKFTextParser())

    def write_to_archive(self):
        self.metainfo_parser.annotation_key = ams.OUT_KEY
        self.archive.data = Simulation(program=Program(name='AMS'))
        self.metainfo_parser.data_object = self.archive.data

        rkf_files = search_files('ams.rkf', os.path.dirname(self.mainfile))
        self.parser = self.mainfile_parser
        self.parser.filepath = self.mainfile
        if rkf_files:
            if len(rkf_files) > 1:
                self.logger.warning('Multiple ams.rkf files found.')
            self.parser = self.rkf_parser
            self.parser.filepath = rkf_files[0]
            self.parser.data_object.parse()

        self.parser.convert(self.metainfo_parser)


class AMSParser(MatchingParser):
    """
    Main parse interface to NOMAD.
    """

    archive_writer = AMSArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
