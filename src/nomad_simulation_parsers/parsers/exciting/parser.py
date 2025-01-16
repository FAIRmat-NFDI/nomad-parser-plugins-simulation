import os
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.parsing.file_parser import Parser
from nomad.parsing.file_parser.mapping_parser import (
    MetainfoParser,
    TextParser,
    XMLParser,
)
from nomad_simulations.schema_packages.general import Simulation

import nomad_simulation_parsers.schema_packages.exciting_schema_package  # noqa
from nomad_simulation_parsers.parsers.utils import search_files

from .eigval_reader import EigvalReader
from .info_reader import InfoReader


class InfoParser(TextParser):
    def get_xc_functionals(self, xc_type: int) -> list[dict[str, Any]]:
        xc_functional_map = {
            2: ['LDA_C_PZ', 'LDA_X_PZ'],
            3: ['LDA_C_PW', 'LDA_X_PZ'],
            4: ['LDA_C_XALPHA'],
            5: ['LDA_C_VBH'],
            20: ['GGA_C_PBE', 'GGA_X_PBE'],
            21: ['GGA_C_PBE', 'GGA_X_PBE_R'],
            22: ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
            26: ['GGA_C_PBE', 'GGA_X_WC'],
            30: ['GGA_C_AM05', 'GGA_C_AM05'],
            300: ['GGA_C_BGCP', 'GGA_X_PBE'],
            406: ['HYB_GGA_XC_PBEH'],
            408: ['HYB_GGA_XC_HSE03'],
        }
        return [dict(libxc=name) for name in xc_functional_map.get(xc_type, [])]

    def get_forces(self, source: dict[str, Any]) -> dict[str, Any]:
        return dict(
            forces=source.get('forces'),
            n_points=len(source.get('forces', [])),
            rank=[3],
        )

    def get_configurations(self, root: dict[str, Any]) -> list[dict[str, Any]]:
        configurations = [
            root[key] for key in ['groundstate', 'hybrid'] if root.get(key)
        ]
        optimization = root.get('structure_optimization')
        if optimization:
            configurations.extend(optimization.get('optimization_step', []))
            configurations.append(optimization)
        return configurations


class InputXMLParser(XMLParser):
    def get_xc_functionals(self, xc_funcs: dict[str, str]) -> list[dict[str, str]]:
        return [dict(libxc=val, type=key) for key, val in xc_funcs.items()]


class EigvalParser(TextParser):
    def get_eigenvalues(self, source: dict[str, Any]):
        eigs_occs = source.get('eigenvalues_occupancies')
        eigs = np.array([v.get('eigenvalues') for v in eigs_occs])
        occs = np.array([v.get('occupancies') for v in eigs_occs])

        return [
            dict(eigenvalues=eigs[:, spin, :], occupancies=occs[:, spin, :])
            for spin in range(len(eigs[0]))
        ]


class ExcitingParser(Parser):
    def parse(
        self, mainfile: str, archive: 'EntryArchive', logger: 'BoundLogger'
    ) -> None:
        maindir = os.path.dirname(mainfile)
        mainbase = os.path.basename(mainfile)

        info_parser = InfoParser(text_parser=InfoReader())
        info_parser.filepath = mainfile

        data_parser = MetainfoParser(data_object=Simulation())
        data_parser.annotation_key = 'info'

        info_parser.convert(data_parser)

        input_xml_files = (
            search_files('input.xml', maindir, mainbase)
            if not archive.m_xpath('data.model_method[0].xc_functionals')
            else []
        )
        if input_xml_files:
            input_xml_parser = InputXMLParser(filepath=input_xml_files[0])
            data_parser.annotation_key = 'input_xml'
            input_xml_parser.convert(data_parser)

        eigval_files = search_files('EIGVAL.OUT', maindir, mainbase)
        if eigval_files:
            eigval_parser = EigvalParser(
                filepath=eigval_files[0], text_parser=EigvalReader()
            )
            data_parser.annotation_key = 'eigval'
            eigval_parser.convert(data_parser, update_mode='merge@-1')
            self.eigval_parser = eigval_parser

        archive.data = data_parser.data_object

        self.info_parser = info_parser
        # close parsers
        # info_parser.close()
        # input_xml_parser.close()
        # eigval_parser.close()
        # data_parser.close()
