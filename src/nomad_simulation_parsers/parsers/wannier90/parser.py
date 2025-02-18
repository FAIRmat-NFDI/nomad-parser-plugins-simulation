import os
import re
from collections.abc import Iterable
from importlib import reload
from typing import Any

import numpy as np
from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import TaskReference
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.parsing.file_parser.text_parser import DataTextParser
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.workflow import SinglePoint
from nomad_simulations.schema_packages.workflow.beyond_dft import DFTTBWorkflow
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import wannier90

from .file_parsers import HrParser, WInParser, WOutParser


configuration = config.get_plugin_entry_point(
    'nomad_simulation_parsers.parsers:wannier90_parser_entry_point'
)


class WHrTextParser(TextParser):
    def get_hoppings(self, source: dict[str, Any], **kwargs) -> dict[str, Any]:
        degeneracy_factors = source.get('degeneracy_factors')[2:]
        full_hoppings = source.get('hoppings', [])
        n_wigner_seitz_points = source.get('degeneracy_factors')[1]
        n_orbitals = source.get('n_orbitals')

        hops = np.reshape(
            full_hoppings,
            (n_wigner_seitz_points, n_orbitals, n_orbitals, 7),
        )

        # storing the crystal field splitting values
        ws0 = int((n_wigner_seitz_points - 1) / 2)
        crystal_fields = [
            hops[ws0, i, i, 5] for i in range(n_orbitals)
        ]  # only real elements

        # delete repeated points for different orbitals
        ws_points = hops[:, :, :, :3]
        ws_points = np.unique(ws_points.reshape(-1, 3), axis=0)

        # passing hoppings
        hoppings = hops[:, :, :, -2] + 1j * hops[:, :, :, -1]
        result = dict(
            degeneracy_factors=degeneracy_factors,
            hoppings=hoppings,
            crystal_fields=crystal_fields,
        )
        if kwargs.get('ws'):
            result.update(dict(ws_points=ws_points, n_ws_points=n_wigner_seitz_points))

        return result


class WDosTextParser(TextParser):
    def get_dos(self, source: np.ndarray) -> dict[str, Any]:
        data = np.transpose(source)
        return dict(energies=data[0], value=data[1])


class WBandTextParser(TextParser):
    def get_data(self, data: np.ndarray) -> np.ndarray:
        return np.transpose(data)[1:].transpose()


class WOutTextParser(TextParser):
    def get_lattice_vectors(self, vectors: list[Any]) -> np.ndarray:
        return np.vstack(vectors[-3:])

    def get_pbc(self, vectors: list[Any]) -> list[bool]:
        return [vectors is not None] * 3

    def is_maximally_localized(self, niter: int, default=0) -> bool:
        return (niter or default) > 1

    def get_kpoints(self, points: np.ndarray) -> np.ndarray:
        return np.complex128(points[::2])

    def get_k_line_path(self, k_line_path: dict[str, Any]):
        high_symm_names = k_line_path.get('high_symm_name')
        high_symm_values = [
            np.reshape(val, (2, 3)) for val in k_line_path.get('high_symm_value')
        ]
        # Start with the first element of the first pair
        names = [high_symm_names[0][0]]
        values = [high_symm_values[0][0]]
        for i, pair in enumerate(high_symm_names):
            # Add the second element if it's not the last one in the list
            if pair[1] != names[-1]:
                names.append(pair[1])
                values.append(high_symm_values[i][1])
        return dict(names=names, values=values)


class WInTextParser(TextParser):
    # TODO these should be defined in common utils
    _l_symbols = ['s', 'p', 'd', 'f']
    _m_symbols = [
        None,
        'x',
        'y',
        'z',
        'z^2',
        'xz',
        'yz',
        'x^2-y^2',
        'xy',
        'z^3',
        'xz^2',
        'yz^2',
        'z(x^2-y^2)',
        'xyz',
        'x(x^2-3y^2)',
        'y(3x^2-y^2)',
    ]
    _wannier_symbols = [
        's',
        'px',
        'py',
        'pz',
        'dz2',
        'dxz',
        'dyz',
        'dx2-y2',
        'dxy',
        'fz3',
        'fxz2',
        'fyz2',
        'fz(x2-y2)',
        'fxyz',
        'fx(x2-3y2)',
        'fy(3x2-y2)',
    ]

    def get_projections(self, source: list[Any]) -> list[dict[str, Any]]:
        return [dict(projection=val) for val in source]

    def get_branch_label_indices(
        self,
        atom: Any,
        positions: list[np.ndarray],
        labels: list[str],
        lattice_vectors: list[np.ndarray],
    ) -> Any:
        symbols, indices = [], []
        if atom is None:
            return None

        elif isinstance(atom, int):
            indices = [atom]

        elif match := re.match(r'([cf])=(.+?),(.+?),(.+)', atom):
            coord = match.groups()[0]
            position = np.array(match.groups()[1:4], float)
            if coord.lower() == 'f':
                position = np.dot(position, lattice_vectors)
            for n, pos in enumerate(positions):
                if np.allclose(
                    position, pos, configuration.equal_cell_positions_tolerance
                ):
                    indices.append(n)
                    symbols.append(labels[n])

        elif isinstance(atom, str):
            indices = [n for n, label in enumerate(labels) if label == atom]
            symbols = [atom]

        return dict(label=''.join(symbols), indices=indices)

    def get_orbitals_state(self, orbital: Any) -> list[dict[str, Any]]:
        if orbital is None:
            return None

        states = []
        orbitals = re.findall(r'l=([\d+])(?:,mr=([\d])+=)?', orbital)
        for orb in orbitals:
            nl = int(orb[0])
            states.append(dict(l=self._l_symbols[nl]))
            if orb[1]:
                nm = sum([len(range(-n, n + 1)) for n in range(nl)]) + int(orb[1])
                states[-1]['m'] = self._m_symbols[nm]
        if not orbitals:
            for orb in orbital.split(';'):
                try:
                    norb = self._wannier_symbols.index(orb)
                except Exception:
                    continue
                # calculate l,m from norb
                nl = 0
                nm = 0
                while True:
                    m_offset = [nm + nq for nq in range(len(range(-nl, nl + 1)))]
                    if norb in m_offset:
                        nm = m_offset.index(norb)
                        break
                    nl += 1
                    nm += len(m_offset)
                states.append(dict(l=self._l_symbols[nl], m=self._m_symbols[nm]))
        return states


class WannierArchiveWriter(ArchiveWriter):
    def parse_workflow(self) -> None:
        """
        Write to archive workflow2 section.
        """
        workflow = SinglePoint()
        workflow.normalize(archive=self.archive, logger=self.logger)
        self.archive.workflow2 = workflow

        # dft+tb workflow
        if self.child_archives:
            dft_file = list(self.child_archives.keys())[0]
            match = re.match(r'.+?/raw/(.+)', dft_file)
            if not match:
                self.logger.error('DFT calculation not found.')
                return
            dft_archive: EntryArchive = self.archive.m_context.resolve_archive(
                f'../upload/archive/mainfile/{match.group(1)}'
            )

            workflow_archive = self.child_archives[dft_file]
            workflow_archive.workflow2 = DFTTBWorkflow(
                tasks=[
                    TaskReference(task=dft_archive.workflow2),
                    TaskReference(task=self.archive.workflow2),
                ]
            )

    def parse_input(self) -> None:
        """
        Parse wannier input files.
        """
        win_files = search_files(
            pattern='*.win', basedir=self.basedir, re_pattern=self.basename
        )
        if not win_files:
            return
        if len(win_files) > 1:
            self.logger.warning('Multiple `*.win` files found, parsing only first.')

        win_parser = WInTextParser(text_parser=WInParser())
        win_parser.filepath = win_files[0]
        # need data from out
        for key in ['structure', 'lattice_vectors']:
            win_parser.data[key] = self.wout_parser.data.get(key)
        self.data_parser.annotation_key = 'win'
        self.data_parser.data_object = self.archive.data
        win_parser.convert(self.data_parser)
        win_parser.close()

    def parse_hr(self) -> None:
        """
        Parse wannier hr files.
        """
        whr_parser = WHrTextParser(text_parser=HrParser())
        hr_files = search_files(
            pattern='*hr.dat', basedir=self.basedir, re_pattern=self.basename
        )
        if len(hr_files) > 1:
            self.logger.info('Multiple `*hr.dat` files found.')
        for hr_file in hr_files:
            whr_parser.filepath = hr_file
            # need data from out
            whr_parser.data['n_orbitals'] = self.wout_parser.data.get('Nwannier')
            self.data_parser.annotation_key = 'whr'
            self.data_parser.data_object = self.archive.data
            whr_parser.convert(self.data_parser)
        whr_parser.close()

    def parse_dos(self) -> None:
        """
        Parse dos files.
        """
        wdos_parser = WDosTextParser(text_parser=DataTextParser())
        dos_files = search_files(
            pattern='*dos.dat', basedir=self.basedir, re_pattern=self.basename
        )
        if len(dos_files) > 1:
            self.logger.info('Multiple `*dos.dat` files found.')
        for dos_file in dos_files:
            wdos_parser.filepath = dos_file
            wdos_parser.data_object.parse('data')
            self.data_parser.annotation_key = 'dos'
            self.data_parser.data_object = self.archive.data
            wdos_parser.convert(self.data_parser)
        wdos_parser.close()

    def parse_band(self) -> None:
        """
        Parse band files.
        """
        wband_parser = WBandTextParser(text_parser=DataTextParser())
        # parse band files
        band_files = search_files(
            pattern='*band.dat', basedir=self.basedir, re_pattern=self.basename
        )
        for band_file in band_files:
            wband_parser.filepath = band_file
            wband_parser.data_object.parse('data')
            self.data_parser.annotation_key = 'band'
            self.data_parser.data_object = self.archive.data
            wband_parser.convert(self.data_parser)
        wband_parser.close()

    def write_to_archive(self) -> None:
        # reload the schema annotations
        reload(wannier90)

        self.basename = os.path.basename(self.mainfile)
        self.basedir = os.path.dirname(self.mainfile)
        # define mapping parser interface to OutParser
        self.wout_parser = WOutTextParser(text_parser=WOutParser())
        self.wout_parser.filepath = self.mainfile

        # construct metainfo parser
        data = Simulation()
        self.data_parser = MetainfoParser()
        self.data_parser.annotation_key = 'wout'
        self.data_parser.data_object = data
        self.archive.data = data

        self.wout_parser.convert(self.data_parser)

        # parse input file
        if data.model_system:
            self.parse_input()

        # parse hr files
        self.parse_hr()

        # parse dos files
        self.parse_dos()

        # parse band files
        self.parse_band()

        self.parse_workflow()

        # close parser contexts
        self.wout_parser.close()
        self.data_parser.close()


class Wannier90Parser(MatchingParser):
    archive_writer = WannierArchiveWriter()
    parser_level = 1

    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable[str]:
        """
        Returns DFT file if present.
        """
        is_mainfile = super().is_mainfile(filename, mime, buffer, decoded_buffer, compression)
        if is_mainfile:
            basedir = os.path.dirname(filename)
            basename = os.path.basename(filename)
            # get the initial DFT calculation
            for f in ['vasprun.xml', 'OUTCAR']:
                files = search_files(
                    pattern=f'*{f}',
                    basedir=basedir,
                    re_pattern=basename,
                )
                if files:
                    self.creates_children = True
                    # return only the first match, this is the key to the
                    # workflow archive
                    return [files[0]]
        return is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = {},
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
