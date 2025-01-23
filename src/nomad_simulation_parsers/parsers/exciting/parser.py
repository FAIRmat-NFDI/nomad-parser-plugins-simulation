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
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.utils.general import (
    remove_mapping_annotations,
    search_files,
)

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

    def get_atoms(self, source: dict[str, Any]) -> dict[str, Any]:
        positions = source.get('positions')
        initial = self.data.get('initialization', {})
        lattice_vectors = initial.get('lattice_vectors')
        if positions is not None and source.get('positions_format') == 'lattice':
            positions = np.dot(positions, lattice_vectors.magnitude)
        if positions is None:
            positions = []
            for species in initial.get('species', []):
                positions_specie = species.get('positions')
                if species.get('positions_format') == 'lattice':
                    positions_specie = np.dot(positions_specie, lattice_vectors)
                positions.extend(positions_specie)
        atoms = []
        exclude = ['positions', 'positions_format', 'radial_points']
        for species in initial.get('species', []):
            atom = {k: v for k, v in species.items() if k not in exclude}
            atoms.extend([atom] * len(species.get('positions', [])))
        if not atoms:
            atoms = [dict(symbol=s) for s in source.get('symbols')]
        return dict(positions=np.array(positions, dtype=float), atoms=atoms)


class InputXMLParser(XMLParser):
    def get_xc_functionals(self, xc_funcs: dict[str, str]) -> list[dict[str, str]]:
        return [dict(libxc=val, type=key) for key, val in xc_funcs.items()]


class BandstructureXMLParser(XMLParser):
    n_spin = 1

    def get_bandstructures(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        # TODO determine format for spin pol case
        energies = [
            p['@eval'] for b in source['bandstructure']['band'] for p in b['point']
        ]
        n_spin = source.get('n_spin', self.n_spin)
        n_band = len(source['bandstructure']['band']) // n_spin
        n_kpoints = len(source['bandstructure']['band'][0]['point'])
        energies = np.array(energies, dtype=float).reshape((n_spin, n_band, n_kpoints))
        return [
            dict(energies=e.T * ureg.hartree, n_states=n_band, n_kpoints=n_kpoints)
            for e in energies
        ]

    def reshape_coords(self, source: list[str]) -> np.ndarray:
        return np.array([v.split() for v in source], dtype=float)


class DosXMLParser(XMLParser):
    def to_float(self, source: list[str]) -> np.ndarray:
        return np.array(source, dtype=float)

    def get_dos(self, source: list[dict[str, Any]]) -> dict[str, Any]:
        return dict(
            dos=np.array([p['@dos'] for p in source.get('point', [])], dtype=float),
            energy=np.array([p['@e'] for p in source.get('point', [])], dtype=float),
        )


class EigvalParser(TextParser):
    def get_eigenvalues(self, source: dict[str, Any]):
        eigs_occs = source.get('eigenvalues_occupancies')
        eigs = np.array([v.get('eigenvalues') for v in eigs_occs])
        occs = np.array([v.get('occupancies') for v in eigs_occs])

        return [
            dict(
                eigenvalues=eigs[:, spin, :],
                occupancies=occs[:, spin, :],
                # n_states printed on file is actual n of states * n spin channels
                n_states=len(eigs[0][spin]),
            )
            for spin in range(len(eigs[0]))
        ]


class ExcitingParser(Parser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = {},
        **kwargs,
    ) -> None:
        from nomad_simulation_parsers.schema_packages import exciting

        maindir = os.path.dirname(mainfile)
        mainbase = os.path.basename(mainfile)

        # mainfile INFO.OUT parser
        info_parser = InfoParser(text_parser=InfoReader())
        info_parser.filepath = mainfile

        data_parser = MetainfoParser(data_object=Simulation())
        data_parser.annotation_key = 'info'

        info_parser.convert(data_parser)

        # read xc functionals from input.xml
        input_xml_files = (
            search_files('input.xml', maindir, mainbase)
            if not archive.m_xpath('data.model_method[0].xc_functionals')
            else []
        )
        if input_xml_files:
            input_xml_parser = InputXMLParser(filepath=input_xml_files[0])
            data_parser.annotation_key = 'input_xml'
            input_xml_parser.convert(data_parser)
            input_xml_parser.close()

        # eigenvalues from eigval.out
        eigval_files = search_files('EIGVAL.OUT', maindir, mainbase)
        if eigval_files:
            eigval_parser = EigvalParser(
                filepath=eigval_files[0], text_parser=EigvalReader()
            )
            data_parser.annotation_key = 'eigval'
            eigval_parser.convert(data_parser, update_mode='merge@-1')
            eigval_parser.close()

        # bandstructure from bandstructure.xml
        bandstructure_files = search_files('bandstructure.xml', maindir, mainbase)
        if bandstructure_files:
            bandstructure_parser = BandstructureXMLParser(
                filepath=bandstructure_files[0]
            )
            # TODO set n_spin from info
            data_parser.annotation_key = 'bandstructure_xml'
            bandstructure_parser.convert(data_parser, update_mode='merge@-1')
            bandstructure_parser.close()

        # dos from dos.xml
        dos_files = search_files('dos.xml', maindir, mainbase)
        if dos_files:
            dos_parser = DosXMLParser(filepath=dos_files[0])
            data_parser.annotation_key = 'dos_xml'
            dos_parser.convert(data_parser, update_mode='merge@-1')
            dos_parser.close()

        archive.data = data_parser.data_object

        # close parsers
        info_parser.close()
        data_parser.close()

        # remove annotations
        remove_mapping_annotations(exciting.general.Simulation.m_def)
