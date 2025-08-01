import os
import re
from collections.abc import Iterable
from typing import Any

import numpy as np
from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import TaskReference
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser, TextParser
from nomad_file_parser.text_parser import DataTextParser
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.workflow import SinglePoint
from nomad_simulations.schema_packages.workflow.dmft import DFTTBDMFTWorkflow
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import wannier90

from .file_parsers import HrParser, WInParser, WOutParser

configuration = config.get_plugin_entry_point(
    'nomad_simulation_parsers.parsers:wannier90_parser'
)
LOGGER = get_logger(__name__)


# TODO temporary fix for structlog unable to propagate logger
class Wannier90MetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class WHrTextParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

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
    # aux data set by the archive writer; kept out of the parsed `data` payload
    _energies_origin = None

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_dos(self, source: np.ndarray) -> dict[str, Any]:
        data = np.transpose(source)
        result = dict(energies=data[0], value=data[1])
        if self._energies_origin is not None:
            result['energies_origin'] = self._energies_origin
        return result


class WBandTextParser(TextParser):
    # aux data set by the archive writer; kept out of the parsed `data` payload
    _k_path = None
    _highest_occupied = None

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_data(self, data: np.ndarray) -> np.ndarray:
        return np.transpose(data)[1:].transpose()

    def get_band_structure(self, data: np.ndarray) -> dict[str, Any]:
        transposed = np.transpose(data)
        band_values = transposed[1:].transpose()

        result = dict(value=band_values)
        if self._k_path:
            result['k_path'] = self._k_path
        if self._highest_occupied is not None:
            result['highest_occupied'] = self._highest_occupied

        return result


class WOutTextParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_lattice_vectors(self, vectors: list[Any]) -> np.ndarray:
        return np.vstack(vectors[-3:])

    def get_pbc(self, vectors: list[Any]) -> list[bool]:
        return [vectors is not None] * 3

    def is_maximally_localized(self, niter: int, default=0) -> bool:
        return (niter or default) > 1

    def get_kpoints(self, points: np.ndarray) -> np.ndarray:
        return np.float64(points[::2])

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
    """
    Parser for Wannier90 .win input files.

    Extracts projection specifications, lattice parameters, and atomic positions
    from Wannier90 input files and converts them to NOMAD schema structures.
    """

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

    # Explicit mapping from Wannier90 symbols to (l, ml) quantum numbers
    # Based on Wannier90 User Guide Table 3.2
    # Real spherical harmonics with Wannier90's specific ordering
    _symbol_to_quantum_numbers = {
        # s orbitals (l=0)
        's': (0, 0),
        # p orbitals (l=1)
        'px': (1, -1),
        'py': (1, 0),
        'pz': (1, 1),
        # d orbitals (l=2)
        'dz2': (2, 0),
        'dxz': (2, 1),
        'dyz': (2, -1),
        'dx2-y2': (2, 2),
        'dxy': (2, -2),
        # f orbitals (l=3) - standard cubic harmonic ordering
        'fz3': (3, 0),
        'fxz2': (3, 1),
        'fyz2': (3, -1),
        'fz(x2-y2)': (3, 2),
        'fxyz': (3, -2),
        'fx(x2-3y2)': (3, 3),
        'fy(3x2-y2)': (3, -3),
    }

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def _get_quantum_numbers_from_symbol(self, symbol: str) -> tuple[int, int] | None:
        """
        Get (l, ml) quantum numbers for a Wannier90 orbital symbol.

        Args:
            symbol: Orbital symbol (e.g., 's', 'px', 'dxy', 'dx2-y2')

        Returns:
            tuple: (l_quantum_number, ml_quantum_number) or None if not found

        Note:
            Uses explicit mapping from Wannier90 User Guide Table 3.2.
            The mapping is for real spherical harmonics with Wannier90's
            specific ordering convention.
        """
        return self._symbol_to_quantum_numbers.get(symbol)

    @staticmethod
    def _calculate_ml_from_position(ll: int, position: int) -> int:
        """
        Convert position within l-manifold to ml quantum number.

        Args:
            ll: Orbital angular momentum quantum number
            position: Position within the l-manifold (0-indexed)

        Returns:
            ml quantum number corresponding to the position

        Note:
            Assumes position 0 → ml=-ll, position 1 → ml=-ll+1, etc.
            This convention should be verified against Wannier90 documentation.
        """
        return -ll + position

    def get_projections(self, source: list[Any]) -> list[dict[str, Any]]:
        return [dict(projection=val) for val in source]

    def get_branch_label_indices(
        self,
        atom: Any,
        positions: list[np.ndarray],
        labels: list[str],
        lattice_vectors: list[np.ndarray],
    ) -> Any:
        """
        Parse atom specification and find matching atom indices in the structure.

        Supports multiple formats:
        - Integer: direct atom index
        - "c=x,y,z": Cartesian coordinates
        - "f=x,y,z": Fractional coordinates
        - String: atom label/symbol

        Args:
            atom: Atom specification from Wannier90 projection
            positions: List of atomic positions
            labels: List of atomic labels/symbols
            lattice_vectors: Unit cell lattice vectors

        Returns:
            Dict with 'label' (combined symbol string) and 'indices'
            (list of atom indices)
        """
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

            # Use global config for tolerance (instance-level tolerance would require
            # refactoring the mapping system to pass parser context)
            tolerance = configuration.equal_cell_positions_tolerance

            for n, pos in enumerate(positions):
                if np.allclose(position, pos, tolerance):
                    indices.append(n)
                    symbols.append(labels[n])

        elif isinstance(atom, str):
            indices = [n for n, label in enumerate(labels) if label == atom]
            symbols = [atom]

        return dict(label=''.join(symbols), indices=indices)

    def get_orbitals_state(self, orbital: Any) -> list[dict[str, Any]]:
        """
        Parse orbital projections from Wannier90 input and create
        ElectronicState dictionary specifications.

        Wannier90 supports two projection formats:
        1. l-based: "l=2,mr=1" - specifies l and position within l-manifold
        2. Symbol-based: "dx2-y2" - uses standard orbital names

        Both formats are converted to ElectronicState structures containing
        SphericalSymmetryState quantum numbers (l, ml).

        Args:
            orbital: Orbital specification string from Wannier90 .win file

        Returns:
            List of dicts with structure:
            [{'spin_orbit_state': {'l_quantum_number': int, 'ml_quantum_number': int}}]
            Returns None if orbital is None.

        Note:
            The ml quantum number mapping assumes Wannier90's ordering convention
            where position within an l-manifold maps linearly to ml values starting
            from ml=-l. This should be verified against Wannier90 documentation.
        """
        if orbital is None:
            return None

        states = []

        # Try parsing l-based format: l=2,mr=1
        orbitals = re.findall(r'l=(\d+)(?:,mr=(\d+)=)?', orbital)
        for orb in orbitals:
            nl = int(orb[0])
            state = {'spin_orbit_state': {'l_quantum_number': nl}}
            if orb[1]:
                # mr parameter specifies position within l-manifold
                ml = self._calculate_ml_from_position(nl, int(orb[1]))
                state['spin_orbit_state']['ml_quantum_number'] = ml
            states.append(state)

        # If l-based format not found, try symbol-based format: dx2-y2
        if not orbitals:
            for orb in orbital.split(';'):
                # Look up quantum numbers from explicit mapping
                quantum_numbers = self._get_quantum_numbers_from_symbol(orb)

                if quantum_numbers is None:
                    # Symbol not recognized, skip it
                    self.logger.warning(
                        f"Unknown orbital symbol '{orb}' in projection specification"
                    )
                    continue

                nl, ml = quantum_numbers
                states.append(
                    {
                        'spin_orbit_state': {
                            'l_quantum_number': nl,
                            'ml_quantum_number': ml,
                        }
                    }
                )

        return states


class WannierArchiveWriter(ArchiveWriter):
    """
    Archive writer for Wannier90 calculations.

    Orchestrates parsing of multiple Wannier90 output files (.wout, .win, _hr.dat,
    _band.dat, etc.) and populates the NOMAD archive with simulation data including
    model systems, methods, and output properties.
    """

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
            if (
                dft_archive is not None
                and dft_archive.workflow2
                and self.archive.workflow2
            ):
                workflow_archive.workflow2 = DFTTBDMFTWorkflow(
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
        self.data_parser.annotation_key = wannier90.WIN_KEY
        self.data_parser.data_object = self.archive.data
        win_parser.convert(self.data_parser)

        reference_energy = win_parser.data.get('energy_fermi')
        if reference_energy is not None:
            self.reference_energy = float(reference_energy) * ureg.eV

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
            self.data_parser.annotation_key = wannier90.WHR_KEY
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
            wdos_parser._energies_origin = self.reference_energy
            self.data_parser.annotation_key = wannier90.DOS_KEY
            self.data_parser.data_object = self.archive.data
            wdos_parser.convert(self.data_parser)
        wdos_parser.close()

    def parse_band(self) -> None:
        """
        Parse band files.
        """
        wband_parser = WBandTextParser(text_parser=DataTextParser())
        k_path = self.wout_parser.get_k_line_path(
            self.wout_parser.data.get('k_line_path')
        )
        # parse band files
        band_files = search_files(
            pattern='*band.dat', basedir=self.basedir, re_pattern=self.basename
        )
        for band_file in band_files:
            wband_parser.filepath = band_file
            wband_parser.data_object.parse('data')
            wband_parser._k_path = k_path
            wband_parser._highest_occupied = self.reference_energy
            self.data_parser.annotation_key = wannier90.BAND_KEY
            self.data_parser.data_object = self.archive.data
            wband_parser.convert(self.data_parser)
        wband_parser.close()

    def write_to_archive(self) -> None:
        self.basename = os.path.basename(self.mainfile)
        self.basedir = os.path.dirname(self.mainfile)
        # define mapping parser interface to OutParser
        self.wout_parser = WOutTextParser(text_parser=WOutParser())
        self.wout_parser.filepath = self.mainfile

        # construct metainfo parser
        data = Simulation()
        self.data_parser = Wannier90MetainfoParser()
        self.data_parser.annotation_key = wannier90.WOUT_KEY
        self.data_parser.data_object = data
        self.archive.data = data
        self.reference_energy = None

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
    """
    Parser for Wannier90 output files.

    Parses Wannier90 calculations including maximally localized Wannier functions,
    band structures, density of states, and hopping matrices.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.archive_writer = WannierArchiveWriter()

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
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
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
                    self.level = 1
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

        # run the old parser
        # TODO remove
        from electronicparsers.wannier90.parser import Wannier90Parser  # noqa

        Wannier90Parser().parse(mainfile, archive, logger)