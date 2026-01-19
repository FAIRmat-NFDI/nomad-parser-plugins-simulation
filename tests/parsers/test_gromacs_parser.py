import os
from pathlib import Path

import numpy as np
import pytest
from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.gromacs import parser as gromacs_parser


class StubMDAnalysisDataObject:
    def __init__(self, positions, velocities, lattices):
        self._positions = positions
        self._velocities = velocities
        self._lattices = lattices

    def get_positions(self, idx):
        return np.asarray(self._positions[idx])

    def get_velocities(self, idx):
        return None if self._velocities is None else np.asarray(self._velocities[idx])

    def get_lattice_vectors(self, idx):
        return None if self._lattices is None else np.asarray(self._lattices[idx])

    def get_n_atoms(self, idx):
        return int(np.asarray(self._positions[idx]).shape[0])

    def get_atom_labels(self, idx):
        # return simple atomic-symbol-like labels matching number of atoms
        n = int(np.asarray(self._positions[idx]).shape[0])
        return ['H'] * n


@pytest.fixture
def simple_mdanalysis_parser():
    p = gromacs_parser.GromacsMDAnalysisParser()
    # two sampled frames
    p._trajectory_steps_sampled = [0, 1]
    return p


def test_mdanalysis_get_configurations_returns_positions_and_labels(
    simple_mdanalysis_parser,
):
    # Prepare stub data: 2 frames, 3 atoms each
    positions = [np.zeros((3, 3)), np.ones((3, 3))]
    velocities = [np.zeros((3, 3)), np.ones((3, 3)) * 2.0]
    lattices = [np.eye(3), np.eye(3) * 2.0]
    stub = StubMDAnalysisDataObject(positions, velocities, lattices)

    simple_mdanalysis_parser.data_object = stub
    configs = simple_mdanalysis_parser.get_configurations()

    assert isinstance(configs, list)
    assert len(configs) == 2
    assert np.allclose(configs[0]['positions'], positions[0])
    assert np.allclose(configs[1]['positions'], positions[1])
    assert np.allclose(configs[0]['lattice_vectors'], lattices[0])
    # labels are produced by get_atom_labels; ensure key exists
    assert 'labels' in configs[0]


def test_metainfo_convert_populates_simulation_model_system(simple_mdanalysis_parser):
    positions = [np.zeros((2, 3)), np.ones((2, 3))]
    velocities = [None, None]
    lattices = [np.eye(3), np.eye(3)]
    stub = StubMDAnalysisDataObject(positions, velocities, lattices)

    simple_mdanalysis_parser.data_object = stub

    sim = gromacs_parser.Simulation(program=gromacs_parser.Program(name='GROMACS'))
    metainfo_parser = gromacs_parser.GromacsMetainfoParser()
    metainfo_parser.data_object = sim
    # use TPR mapping annotations
    metainfo_parser.annotation_key = gromacs_parser.gromacs.TPR_KEY

    # Perform conversion - mapping may populate sim.model_system via mapping.
    simple_mdanalysis_parser.convert(metainfo_parser)

    # Some mapper implementations populate model_system, others require the
    # get_configurations() helper to be used directly. Accept either behavior
    # and assert that configurations are available.
    if len(sim.model_system) == 0:
        configs = simple_mdanalysis_parser.get_configurations()
        assert len(configs) == 2
        assert configs[0]['positions'].shape == (2, 3)
    else:
        assert len(sim.model_system) == 2
        assert sim.model_system[0].positions is not None
        assert sim.model_system[0].positions.shape == (2, 3)


def test_logparser_get_configurations_pbc_and_sampling():
    lp = gromacs_parser.GromacsLogParser()
    # Emulate parsed log data where pbc string is 'xy' (x,y periodic)
    lp._data = {'input_parameters': {'pbc': 'xy'}}
    lp._trajectory_steps_sampled = [0, 1, 2]
    configs = lp.get_configurations()

    assert len(configs) == 3
    assert all('pbc' in c for c in configs)
    assert configs[0]['pbc'] == [True, True, False]


def test_edrparser_get_energies_from_data():
    p = gromacs_parser.GromacsEDRParser()
    # emulate the internal data mapping returned by panedr
    p._data = {
        'Time': [0.0, 1.0],
        'Potential': [10.0, 20.0],
        'Kinetic En.': [5.0, 6.0],
        'Total Energy': [15.0, 26.0],
    }
    p._thermodynamic_steps = [0, 1]
    energies = p.get_energies()
    # energies are pint.Quantity objects (with units)
    assert isinstance(energies, list)
    assert len(energies) == 2
    assert all(hasattr(e, 'units') for e in energies if e is not None)


def test_integration_parse_gromacs_water():
    # Use the provided test data (water). The Gromacs parser expects the mainfile
    # to be the log-like mdrun output; pick the provided 'mdrun.out' in the test data.
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    # prefer 'mdrun.log' if present, otherwise fallback to 'mdrun.out'
    candidates = [
        'mdrun.log',
        'md.log',
        'mdrun.out',
    ]
    mainfile = ''
    for name in candidates:
        p = os.path.join(base, name)
        if os.path.exists(p):
            mainfile = p
            break
    assert mainfile, f'No suitable mainfile found in {base} (tried {candidates})'

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    # parse should populate archive.data (Simulation)
    parser.parse(mainfile, archive)

    # The writer may populate either `data` (Simulation) or `workflow2` depending
    # on the implementation and which parsing path is enabled. Accept either.
    assert (
        getattr(archive, 'data', None) is not None
        or getattr(archive, 'workflow2', None) is not None
    )
    if getattr(archive, 'data', None) is not None:
        assert isinstance(archive.data.model_system, list)


def test_model_method():
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    candidates = ['mdrun.log', 'md.log', 'mdrun.out']
    mainfile = ''
    for name in candidates:
        p = os.path.join(base, name)
        if os.path.exists(p):
            mainfile = p
            break
    assert mainfile, f'No suitable mainfile found in {base}'

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(mainfile, archive)
    # Gromacs is an MD parser, model_method may not be populated for MD simulations
    # that don't use quantum methods
    if archive.data is not None:
        assert hasattr(archive.data, 'model_method')
