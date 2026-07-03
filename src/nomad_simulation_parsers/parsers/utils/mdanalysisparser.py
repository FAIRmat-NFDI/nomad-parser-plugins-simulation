#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from __future__ import annotations

import os
from array import array
from collections import namedtuple
from typing import TYPE_CHECKING, Any, NamedTuple

import numpy as np
from scipy import sparse
from scipy.stats import linregress

try:
    import MDAnalysis
    import MDAnalysis.analysis.rdf as MDA_RDF
    from MDAnalysis.topology.guessers import guess_atom_element

    _HAS_MDA = True
except ImportError:
    _HAS_MDA = False
    MDAnalysis = None  # type: ignore

# Import types for type checking only (not at runtime)
if TYPE_CHECKING:
    from MDAnalysis import Universe as MDAUniverse

from nomad.parsing.file_parser import FileParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.utils.molecular_dynamics import (
    BeadGroup,
    shifted_correlation_average,
)

from nomad_simulation_parsers.parsers.utils.constants import MOLE

LOGGER = get_logger(__name__)


def _check_mda_dependency(function_name: str) -> bool:
    """
    Checks if MDAnalysis is available and logs an error if not.

    Args:
        function_name: The name of the function requiring MDAnalysis

    Returns:
        True if MDAnalysis is available, False otherwise
    """
    if not _HAS_MDA:
        LOGGER.error(
            f'{function_name} requires MDAnalysis. '
            'Please install with `pip install nomad-simulation-parsers[md]`.'
        )
        return False
    return True


class MDAnalysisParser(FileParser):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self._args = args
        self._kwargs = kwargs
        self._atomsgroup_info = None
        self._results = None
        self.universe_error = None

    @property
    def auxilliary_files(self):
        return self._args

    @auxilliary_files.setter
    def auxilliary_files(self, value):
        self._file_handler = None
        self._args = [value] if isinstance(value, str) else value

    @property
    def options(self):
        return self._kwargs

    @options.setter
    def options(self, value):
        self._file_handler = None
        self._kwargs = value

    @property
    def universe(self) -> MDAUniverse | None:
        if not _check_mda_dependency('universe'):
            return None
        # Avoid repeatedly retrying and re-logging after a known init failure.
        if self._file_handler is None and self.universe_error is not None:
            return None
        if self._file_handler is None:
            try:
                self._file_handler = MDAnalysis.Universe(
                    self.mainfile, *self.auxilliary_files, **self.options
                )
            except Exception as e:
                self.logger.error('Error creating MDAnalysis universe.', exc_info=e)
                self.universe_error = e
        return self._file_handler

    @property
    def bead_groups(self):
        atoms_moltypes = self.get('atoms_info', {}).get('moltypes', [])
        moltypes = np.unique(atoms_moltypes)
        bead_groups = {}
        compound = 'fragments'
        for moltype in moltypes:
            if hasattr(self.universe.atoms, 'moltypes'):
                ags_by_moltype = self.universe.select_atoms('moltype ' + moltype)
            else:  # this is easier than adding something to the universe
                selection = ' '.join(
                    [str(i) for i in np.where(atoms_moltypes == moltype)[0]]
                )
                selection = f'index {selection}'
                ags_by_moltype = self.universe.select_atoms(selection)
            ags_by_moltype = ags_by_moltype[ags_by_moltype.masses > abs(1e-2)]
            # remove any virtual/massless sites (needed for, e.g., 4-bead water models)
            bead_groups[moltype] = BeadGroup(ags_by_moltype, compound=compound)

        return bead_groups

    def parse(self, quantity_key: str = None, **kwargs):
        if self._results is None:
            self._results: dict[str, Any] = dict()

        if not self.universe:
            return

        atoms = list(self.universe.atoms)

        name_map = {'mass': 'masses'}
        unit_map = {'mass': ureg.amu, 'charge': ureg.elementary_charge}
        self._results['atoms_info'] = dict()
        for key in [
            'name',
            'charge',
            'mass',
            'resid',
            'resname',
            'molnum',
            'moltype',
            'type',
            'segid',
            'element',
        ]:
            try:
                value = [getattr(atom, key) for atom in atoms]
            except Exception:
                value = None
            value = value * unit_map.get(key, 1) if value is not None else value
            self._results['atoms_info'][name_map.get(key, f'{key}s')] = value

        # Fallbacks/substitutions for missing atom info
        substitutions = [
            ('names', lambda: ['CGX'] * self.universe.atoms.n_atoms),
            (
                'moltypes',
                lambda: (
                    self.get_fragtypes()
                    if hasattr(self.universe.atoms, 'fragments')
                    else None
                ),
            ),
            ('molnums', lambda: getattr(self.universe.atoms, 'fragindices', None)),
            ('resnames', lambda: self._results['atoms_info'].get('resids')),
            ('names', lambda: self._results['atoms_info'].get('types')),
            ('elements', lambda: self._results['atoms_info'].get('names')),
        ]
        for key, func in substitutions:
            if self._results['atoms_info'].get(key) is None:
                try:
                    val = func()
                    if val is not None:
                        self._results['atoms_info'][key] = val
                except Exception:
                    continue

        self._results['n_particles'] = self.universe.atoms.n_atoms
        self._results['n_frames'] = len(self.universe.trajectory)

    def get_fragtypes(self):
        # TODO put description otherwise, make private or put under parse method
        """ """
        atoms_fragtypes = np.empty(self.universe.atoms.types.shape, dtype=str)
        ctr_fragtype = 0
        atoms_fragtypes[self.universe.atoms.fragments[0].ix] = ctr_fragtype
        frag_unique_atomtypes = [
            self.universe.atoms.types[self.universe.atoms.fragments[0].ix]
        ]
        ctr_fragtype += 1
        for i_frag in range(1, self.universe.atoms.n_fragments):
            types_i_frag = self.universe.atoms.types[
                self.universe.atoms.fragments[i_frag].ix
            ]
            flag_fragtype_exists = False
            sorted_i = np.sort(types_i_frag)
            for j_frag in range(len(frag_unique_atomtypes) - 1, -1, -1):
                types_j_frag = frag_unique_atomtypes[j_frag]
                if len(types_i_frag) != len(types_j_frag):
                    continue
                elif np.all(sorted_i == np.sort(types_j_frag)):
                    atoms_fragtypes[self.universe.atoms.fragments[i_frag].ix] = j_frag
                    flag_fragtype_exists = True
            if not flag_fragtype_exists:
                atoms_fragtypes[self.universe.atoms.fragments[i_frag].ix] = ctr_fragtype
                frag_unique_atomtypes.append(
                    self.universe.atoms.types[self.universe.atoms.fragments[i_frag].ix]
                )
                ctr_fragtype += 1
        return atoms_fragtypes

    class RDFParams(NamedTuple):
        n_traj_split: int = 10
        n_prune: int = 1
        interval_indices: list[list[int]] | None = None
        n_bins: int = 200
        n_smooth: int = 2
        max_mols: int = 5000

    def calc_molecular_rdf(self, params: RDFParams | None = None):
        """
        Calculates the radial distribution functions between for each unique pair of
        molecule types as a function of their center of mass distance.

        Args:
            params (RDFParams, optional): Parameters for RDF calculation.
        """
        params = params or self.RDFParams()
        if self.universe is None:
            return
        trajectory = self.universe.trajectory[0] if self.universe.trajectory else None
        dimensions = getattr(trajectory, 'dimensions', None) if trajectory else None
        if dimensions is None:
            return

        n_frames = self.universe.trajectory.n_frames
        frames_start, frames_end, n_frames_split, interval_indices = (
            self._get_rdf_intervals(
                n_frames, params.n_traj_split, params.interval_indices
            )
        )

        bead_groups = self.bead_groups
        if not bead_groups:
            return bead_groups
        moltypes = [moltype for moltype in bead_groups.keys()]
        moltypes = self._prune_large_rdf_bead_groups(
            moltypes, bead_groups, params.max_mols
        )

        min_box_dimension = np.min(self.universe.trajectory[0].dimensions[:3])
        max_rdf_dist = min_box_dimension / 2

        rdf_results = self._initialize_rdf_results(params.n_smooth)
        for i, moltype_i in enumerate(moltypes):
            for j, moltype_j in enumerate(moltypes):
                if not self._should_calculate_rdf(i, j, bead_groups, moltype_i):
                    continue
                exclusion_block = (1, 1) if i == j else None
                pair_type = f'{moltype_i}-{moltype_j}'
                rdf_pair_params = self.RDFPairParams(
                    bead_groups=bead_groups,
                    moltype_i=moltype_i,
                    moltype_j=moltype_j,
                    pair_type=pair_type,
                    params=params,
                    frames_start=frames_start,
                    frames_end=frames_end,
                    max_rdf_dist=max_rdf_dist,
                    exclusion_block=exclusion_block,
                )
                rdf_results_tmp = self._calculate_rdf_for_pair(rdf_pair_params)
                for interval_group in interval_indices:
                    self._get_rdf_avg(
                        rdf_results_tmp, rdf_results, interval_group, n_frames_split
                    )
        return rdf_results

    def _get_rdf_intervals(self, n_frames, n_traj_split, interval_indices):
        if n_frames < n_traj_split:
            n_traj_split = 1
            frames_start = np.array([0])
            frames_end = np.array([n_frames])
            n_frames_split = np.array([n_frames])
            interval_indices = [[0]]
        else:
            run_len = int(n_frames / n_traj_split)
            frames_start = np.arange(n_traj_split) * run_len
            frames_end = frames_start + run_len
            frames_end[-1] = n_frames
            n_frames_split = frames_end - frames_start
            if not interval_indices:
                interval_indices = [[i] for i in range(n_traj_split)]
        return frames_start, frames_end, n_frames_split, interval_indices

    def _prune_large_rdf_bead_groups(self, moltypes, bead_groups, max_mols):
        del_list = []
        for i_moltype, moltype in enumerate(moltypes):
            if bead_groups[moltype]._nbeads > max_mols:
                del_list.append(i_moltype)
                self.logger.warning(
                    'The number of molecules exceeds the maximum for calculating '
                    'the rdf. Skipping this molecule type.'
                )
        return np.delete(moltypes, del_list)

    def _initialize_rdf_results(self, n_smooth):
        return {
            'n_smooth': n_smooth,
            'types': [],
            'variables_name': [],
            'bins': [],
            'value': [],
            'frame_start': [],
            'frame_end': [],
        }

    def _should_calculate_rdf(self, i, j, bead_groups, moltype_i):
        if j > i:
            return False
        if i == j and bead_groups[moltype_i].positions.shape[0] == 1:
            return False
        return True

    class RDFPairParams(NamedTuple):
        bead_groups: dict
        moltype_i: Any
        moltype_j: Any
        pair_type: str
        params: Any
        frames_start: np.ndarray
        frames_end: np.ndarray
        max_rdf_dist: float
        exclusion_block: Any

    def _calculate_rdf_for_pair(self, rdf_pair_params: MDAnalysisParser.RDFPairParams):
        bead_groups = rdf_pair_params.bead_groups
        moltype_i = rdf_pair_params.moltype_i
        moltype_j = rdf_pair_params.moltype_j
        pair_type = rdf_pair_params.pair_type
        params = rdf_pair_params.params
        frames_start = rdf_pair_params.frames_start
        frames_end = rdf_pair_params.frames_end
        max_rdf_dist = rdf_pair_params.max_rdf_dist
        exclusion_block = rdf_pair_params.exclusion_block

        rdf_results_tmp = {
            'types': [],
            'variables_name': [],
            'bins': [],
            'value': [],
            'frame_start': [],
            'frame_end': [],
        }
        for i_interval in range(params.n_traj_split):
            rdf_results_tmp['types'].append(pair_type)
            rdf_results_tmp['variables_name'].append(['distance'])
            rdf = MDA_RDF.InterRDF(
                bead_groups[moltype_i],
                bead_groups[moltype_j],
                range=(0, max_rdf_dist),
                exclusion_block=exclusion_block,
                nbins=params.n_bins,
            ).run(frames_start[i_interval], frames_end[i_interval], params.n_prune)
            rdf_results_tmp['frame_start'].append(frames_start[i_interval])
            rdf_results_tmp['frame_end'].append(frames_end[i_interval])
            n_smooth = params.n_smooth
            rdf_results_tmp['bins'].append(
                rdf.results.bins[int(n_smooth / 2) : -int(n_smooth / 2)] * ureg.angstrom
            )
            rdf_results_tmp['value'].append(
                np.convolve(
                    rdf.results.rdf,
                    np.ones((n_smooth,)) / n_smooth,
                    mode='same',
                )[int(n_smooth / 2) : -int(n_smooth / 2)]
            )
        return rdf_results_tmp

    def _get_rdf_avg(
        self, rdf_results_tmp, rdf_results, interval_indices, n_frames_split
    ):
        RDF_SPLIT_WEIGHT_TOL = 1e-6
        split_weights = n_frames_split[np.array(interval_indices)] / np.sum(
            n_frames_split[np.array(interval_indices)]
        )
        assert abs(np.sum(split_weights) - 1.0) < RDF_SPLIT_WEIGHT_TOL
        rdf_values_avg = (
            split_weights[0] * rdf_results_tmp['value'][interval_indices[0]]
        )
        for i_interval, interval in enumerate(interval_indices[1:]):
            assert (
                rdf_results_tmp['types'][interval]
                == rdf_results_tmp['types'][interval - 1]
            )
            assert (
                rdf_results_tmp['variables_name'][interval]
                == rdf_results_tmp['variables_name'][interval - 1]
            )
            assert (
                rdf_results_tmp['bins'][interval]
                == rdf_results_tmp['bins'][interval - 1]
            ).all()
            rdf_values_avg += (
                split_weights[i_interval + 1] * rdf_results_tmp['value'][interval]
            )
        rdf_results['types'].append(rdf_results_tmp['types'][interval_indices[0]])
        rdf_results['variables_name'].append(
            rdf_results_tmp['variables_name'][interval_indices[0]]
        )
        rdf_results['bins'].append(rdf_results_tmp['bins'][interval_indices[0]])
        rdf_results['value'].append(rdf_values_avg)
        rdf_results['frame_start'].append(
            int(rdf_results_tmp['frame_start'][interval_indices[0]])
        )
        rdf_results['frame_end'].append(
            int(rdf_results_tmp['frame_end'][interval_indices[-1]])
        )

    @property
    def with_trajectory(self):
        """
        True if trajectory is present.
        """
        universe = self.universe
        return (
            universe is not None
            and universe.trajectory is not None
            and len(universe.trajectory) > 0
        )

    def get_frame(self, frame_index):
        """
        Returns the frame in the trajectory with index frame_index.
        """
        universe = self.universe
        if universe is None:
            return None
        try:
            return universe.trajectory[frame_index]
        except Exception as e:
            self.logger.warning('Error accessing frame.', exc_info=e)
            return None

    def get_n_atoms(self, frame_index):
        """
        Returns the number of atoms of the frame with index frame_index.
        """
        frame = self.get_frame(frame_index)
        return len(frame) if frame is not None else None

    def get_atom_labels(self, frame_index):
        """
        Returns the number of atoms of the frame with index frame_index.
        """
        # MDAnalysis assumes no change in atom configuration
        return [
            guess_atom_element(name).title()
            for name in self.get('atoms_info', {}).get('names', [])
        ]

    def get_time(self, frame_index):
        """
        Returns the elapsed simulated physical time since the start of the simulation
        for index frame_index.
        """
        frame = self.get_frame(frame_index)
        return frame.time * ureg.picosecond if frame is not None else None

    def get_step(self, frame_index):
        """
        Returns the step of the frame with index frame_index.
        """
        frame = self.get_frame(frame_index)
        if frame is None:
            return None
        dt = frame.dt if frame.dt else getattr(self.universe.trajectory, 'dt', None)
        if not dt:
            return
        return round(frame.time / dt)

    def get_lattice_vectors(self, frame_index):
        """
        Returns the lattice vectors of the frame with index frame_index.
        """
        return self.get_frame_data(frame_index)['lattice_vectors']

    def get_pbc(self, frame_index):
        """
        Returns the lattice periodicity of the frame with index frame_index.
        """
        lattice_vectors = self.get_lattice_vectors(frame_index)
        return [True] * 3 if lattice_vectors is not None else [False] * 3

    def get_frame_data(self, frame_index) -> dict:
        """
        Returns positions, velocities and lattice_vectors for frame_index in a single
        trajectory seek. Callers that need all three quantities should prefer this over
        calling get_positions/get_velocities/get_lattice_vectors individually.
        """
        frame = self.get_frame(frame_index)
        if frame is None:
            return dict(positions=None, velocities=None, lattice_vectors=None)
        lv = frame.triclinic_dimensions
        return dict(
            positions=frame.positions * ureg.angstrom if frame.has_positions else None,
            velocities=(
                frame.velocities * ureg.angstrom / ureg.ps
                if frame.has_velocities
                else None
            ),
            lattice_vectors=lv * ureg.angstrom if lv is not None else None,
        )

    def get_positions(self, frame_index):
        """
        Returns the positions of the atoms of the frame with index frame_index.
        """
        return self.get_frame_data(frame_index)['positions']

    def get_velocities(self, frame_index):
        """
        Returns the velocities of the atoms of the frame with index frame_index.
        """
        return self.get_frame_data(frame_index)['velocities']

    def get_forces(self, frame_index):
        """
        Returns the forces on the atoms of the frame with index frame_index.
        """
        frame = self.get_frame(frame_index)
        if frame is None:
            return None
        return (
            frame.forces * ureg.kJ / (MOLE * ureg.angstrom)
            if frame.has_forces
            else None
        )

    def get_interactions(self):
        interactions = self.get('interactions', None)
        if interactions is not None:
            return interactions

        if self._results is None:
            self._results = {}

        interaction_types = ['angles', 'bonds', 'dihedrals', 'impropers']
        interactions = []
        for interaction_type in interaction_types:
            try:
                interaction = getattr(self.universe, interaction_type)
            except Exception:
                continue

            for inter in interaction:
                atom_labels = None
                try:
                    atom_labels = [
                        self.universe.atoms[ind].type for ind in inter.indices
                    ]
                except Exception:
                    self.logger.warning('Could not assign atom labels to interactions.')
                interactions.append(
                    dict(
                        atom_labels=atom_labels,
                        # parameters=float(inter.value()),
                        ## This is not the parameter but rather the value of the
                        ## interaction order parameter for a single frame
                        # TODO implement functions to get parameters for
                        # TODO individual parsers
                        atom_indices=inter.indices,
                        type=inter.btype,
                    )
                )

        self._results['interactions'] = interactions

        return interactions

    def __calc_diffusion_constant(
        self, times: np.ndarray, values: np.ndarray, dim: int = 3
    ):
        """
        Determines the diffusion constant from a fit of the mean squared displacement
        vs. time according to the Einstein relation.
        """
        linear_model = linregress(times, values)
        slope = linear_model.slope
        error = linear_model.rvalue
        return slope * 1 / (2 * dim), error

    MIN_FRAMES_FOR_MSD = 50

    def calc_molecular_mean_squared_displacements(self, max_mols=5000):
        """
        Calculates the mean squared displacement for the center of mass of each
        molecule type.

        max_mols: the maximum number of molecules per bead group for calculating
          the msd, for efficiency purposes. 1M is arbitrary, 50k was tested and is
          very fast and does not seem to have any memory issues.
        """
        if self.universe is None:
            return
        trajectory = self.universe.trajectory[0] if self.universe.trajectory else None
        dimensions = getattr(trajectory, 'dimensions', None) if trajectory else None
        if dimensions is None:
            return

        n_frames = self.universe.trajectory.n_frames
        if n_frames < self.MIN_FRAMES_FOR_MSD:
            self.logger.warning(
                'Not enough frames to calculate molecular'
                ' mean squared displacements, skipping.',
            )
            return

        dt = getattr(self.universe.trajectory, 'dt')
        if dt is None:
            return
        times = np.arange(n_frames) * dt

        bead_groups = self.bead_groups
        if bead_groups == {}:
            return bead_groups

        moltypes = [moltype for moltype in bead_groups.keys()]
        moltypes, bead_groups = self._prune_large_bead_groups(
            moltypes, bead_groups, max_mols
        )

        msd_results = self._calculate_msd_for_moltypes(moltypes, bead_groups, times)
        return msd_results

    def _prune_large_bead_groups(self, moltypes, bead_groups, max_mols):
        del_list = []
        for i_moltype, moltype in enumerate(moltypes):
            if bead_groups[moltype]._nbeads > max_mols:
                try:
                    moltype_indices = np.array(
                        [atom._ix for atom in bead_groups[moltype]._atoms]
                    )
                    molnums = self.universe.atoms.molnums[moltype_indices]
                    molnum_types = np.unique(molnums)
                    molnum_types_rnd = np.sort(
                        np.random.choice(molnum_types, size=max_mols)
                    )
                    atom_indices_rnd = moltype_indices[
                        np.concatenate(
                            [
                                np.where(molnums == molnum)[0]
                                for molnum in molnum_types_rnd
                            ]
                        )
                    ]
                    selection = ' '.join([str(i) for i in atom_indices_rnd])
                    selection = f'index {selection}'
                    ags_moltype_rnd = self.universe.select_atoms(selection)
                    bead_groups[moltype] = BeadGroup(
                        ags_moltype_rnd, compound='fragments'
                    )
                    self.logger.warning(
                        'Maximum number of molecules for calculating the msd has been '
                        'reached. Will make a random selection for calculation.'
                    )
                except Exception:
                    self.logger.warning(
                        'Tried to select random molecules for large group when '
                        'calculating msd, but something went wrong. '
                        'Skipping this molecule type.'
                    )
                del_list.append(i_moltype)
        moltypes = np.delete(moltypes, del_list)
        return moltypes, bead_groups

    def _mean_squared_displacement(self, start: np.ndarray, current: np.ndarray):
        """
        Calculates mean square displacement between current and initial (start)
        coordinates.
        """
        vec = start - current
        return (vec**2).sum(axis=1).mean()

    def _calculate_msd_for_moltypes(self, moltypes, bead_groups, times):
        msd_results = {}
        msd_results['value'] = []
        msd_results['times'] = []
        msd_results['diffusion_constant'] = []
        msd_results['error_diffusion_constant'] = []
        for moltype in moltypes:
            positions = self.get_nojump_positions(bead_groups[moltype])
            results = shifted_correlation_average(
                self._mean_squared_displacement, times, positions
            )
            msd_results['value'].append(results[1])
            msd_results['times'].append(results[0])
            diffusion_constant, error = self.__calc_diffusion_constant(*results)
            msd_results['diffusion_constant'].append(diffusion_constant)
            msd_results['error_diffusion_constant'].append(error)

        msd_results['types'] = moltypes
        msd_results['times'] = np.array(msd_results['times']) * ureg.picosecond
        msd_results['value'] = np.array(msd_results['value']) * ureg.angstrom**2
        msd_results['diffusion_constant'] = (
            np.array(msd_results['diffusion_constant'])
            * ureg.angstrom**2
            / ureg.picosecond
        )
        msd_results['error_diffusion_constant'] = np.array(
            msd_results['error_diffusion_constant']
        )
        return msd_results

    def parse_jumps(self, selection):
        __ = self.universe.trajectory[0]
        prev = np.array(selection.positions)
        box = self.universe.trajectory[0].dimensions[:3]
        sparse_data = namedtuple('SparseData', ['data', 'row', 'col'])
        jump_data = (
            sparse_data(data=array('b'), row=array('l'), col=array('l')),
            sparse_data(data=array('b'), row=array('l'), col=array('l')),
            sparse_data(data=array('b'), row=array('l'), col=array('l')),
        )

        for i_frame, _ in enumerate(self.universe.trajectory[1:]):
            curr = np.array(selection.positions)
            delta = ((curr - prev) / box).round().astype(np.int8)
            prev = np.array(curr)
            for d in range(3):
                (col,) = np.where(delta[:, d] != 0)
                jump_data[d].col.extend(col)
                jump_data[d].row.extend([i_frame] * len(col))
                jump_data[d].data.extend(delta[col, d])

        return jump_data

    def generate_nojump_matrices(self, selection):
        jump_data = self.parse_jumps(selection)
        N = len(self.universe.trajectory)
        M = selection.positions.shape[0]

        nojump_matrices = tuple(
            sparse.csr_matrix((np.array(m.data), (m.row, m.col)), shape=(N, M))
            for m in jump_data
        )
        return nojump_matrices

    def get_nojump_positions(self, selection):
        nojump_matrices = self.generate_nojump_matrices(selection)
        box = self.universe.trajectory[0].dimensions[:3]

        nojump_positions = []
        for i_frame, __ in enumerate(self.universe.trajectory):
            delta = (
                np.array(
                    np.vstack([m[:i_frame, :].sum(axis=0) for m in nojump_matrices]).T
                )
                * box
            )
            nojump_positions.append(selection.positions - delta)

        return np.array(nojump_positions)

    def clean(self):
        for name in os.listdir(self.maindir):
            if name.startswith('.') and (
                name.endswith('.lock') or name.endswith('.npz')
            ):
                os.remove(os.path.join(self.maindir, name))
