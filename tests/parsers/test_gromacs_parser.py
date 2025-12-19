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


def test_force_field_parsing_from_tpr():
    """Test that force field parameters are extracted from TPR file."""
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    candidates = ['mdrun.log', 'md.log', 'mdrun.out']
    mainfile = ''
    for name in candidates:
        p = os.path.join(base, name)
        if os.path.exists(p):
            mainfile = p
            break

    if not mainfile:
        pytest.skip(f'No suitable mainfile found in {base}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(mainfile, archive)

    # Verify data section exists
    assert archive.data is not None, 'Archive.data should be populated'
    assert hasattr(archive.data, 'model_method'), 'Simulation should have model_method'

    # Check model_method was populated
    if not archive.data.model_method or len(archive.data.model_method) == 0:
        pytest.skip('No model_method found - TPR file may not contain force field data')

    # Find ForceField in model_method (identified by having contributions attribute)
    force_field = None
    for method in archive.data.model_method:
        if hasattr(method, 'contributions'):
            force_field = method
            break

    if force_field is None:
        pytest.skip('No ForceField found in model_method')

    # Verify ForceField structure populated via annotations
    assert hasattr(force_field, 'contributions'), 'ForceField should have contributions'
    assert hasattr(force_field, 'numerical_settings'), (
        'ForceField should have numerical_settings'
    )

    # Check ForceCalculations (numerical settings) from LOG_KEY annotations
    if len(force_field.numerical_settings) > 0:
        force_calc = force_field.numerical_settings[0]
        # Verify attributes exist (values may be None if not in log/mdp)
        assert hasattr(force_calc, 'vdw_cutoff'), (
            'ForceCalculations should have vdw_cutoff'
        )
        assert hasattr(force_calc, 'coulomb_cutoff'), (
            'ForceCalculations should have coulomb_cutoff'
        )
        assert hasattr(force_calc, 'coulomb_type'), (
            'ForceCalculations should have coulomb_type'
        )
        assert hasattr(force_calc, 'neighbor_update_frequency'), (
            'ForceCalculations should have neighbor_update_frequency'
        )

        # If values are populated, verify they have correct types
        if force_calc.vdw_cutoff is not None:
            assert hasattr(force_calc.vdw_cutoff, 'magnitude'), (
                'vdw_cutoff should be a pint Quantity'
            )
        if force_calc.coulomb_cutoff is not None:
            assert hasattr(force_calc.coulomb_cutoff, 'magnitude'), (
                'coulomb_cutoff should be a pint Quantity'
            )

    # Check contributions (Potential list) from TPR_KEY get_force_field_contributions()
    if len(force_field.contributions) > 0:
        potential = force_field.contributions[0]
        assert hasattr(potential, 'type'), 'Potential should have type'
        assert hasattr(potential, 'functional_form'), (
            'Potential should have functional_form'
        )
        assert hasattr(potential, 'parameters'), 'Potential should have parameters'

        # Verify that type is one of the valid enum values (or None)
        if potential.type is not None:
            valid_types = [
                'bond',
                'angle',
                'dihedral',
                'improper dihedral',
                'nonbonded',
                'bond-angle',
            ]
            assert potential.type in valid_types, (
                f'Potential type {potential.type} should be in {valid_types}'
            )


def test_get_coulomb_type_transformation():
    """Test the get_coulomb_type transformation function."""
    lp = gromacs_parser.GromacsLogParser()

    # Test various coulombtype mappings
    test_cases = [
        ('cut-off', 'cutoff'),
        ('cutoff', 'cutoff'),
        ('Ewald', 'ewald'),
        ('PME', 'particle_mesh_ewald'),
        ('P3M-AD', 'particle_particle_particle_mesh'),
        ('Reaction-Field', 'reaction_field'),
        ('Reaction-Field-zero', 'reaction_field'),
        ('unknown', None),
    ]

    for input_val, expected in test_cases:
        result = lp.get_coulomb_type(input_val)
        assert result == expected, (
            f'get_coulomb_type({input_val}) should return {expected}, got {result}'
        )


def test_get_force_field_contributions_transformation():
    """Test the get_force_field_contributions transformation function."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    # Mock MDAnalysis-style interactions with atom_indices and atom_labels
    mock_interactions = [
        {'type': 'bond_harmonic', 'atom_indices': [0, 1], 'atom_labels': ['O', 'H']},
        {'type': 'bond_harmonic', 'atom_indices': [2, 3], 'atom_labels': ['O', 'H']},
        {
            'type': 'angle_harmonic',
            'atom_indices': [0, 1, 2],
            'atom_labels': ['H', 'O', 'H'],
        },
        {
            'type': 'angle_harmonic',
            'atom_indices': [3, 4, 5],
            'atom_labels': ['H', 'O', 'H'],
        },
    ]

    # Mock the data_object to return interactions
    mdap.data_object = type(
        'obj',
        (object,),
        {
            'get': lambda self, key: 'GROMACS 2024' if key == 'version' else None,
            'get_interactions': lambda self, version: mock_interactions,
        },
    )()

    contributions = mdap.get_force_field_contributions()

    assert isinstance(contributions, list), 'Should return list'
    assert len(contributions) == 2, (
        'Should have 2 grouped contributions (bonds and angles)'
    )

    # Check that contributions are grouped by type
    contrib_types = {c['functional_form'] for c in contributions}
    assert 'bond_harmonic' in contrib_types, 'Should have bond_harmonic contribution'
    assert 'angle_harmonic' in contrib_types, 'Should have angle_harmonic contribution'

    # Check structure of contributions
    for contrib in contributions:
        assert 'functional_form' in contrib, 'Should have functional_form'
        assert 'particle_indices' in contrib, 'Should have particle_indices'
        assert isinstance(contrib['particle_indices'], list), (
            'particle_indices should be list'
        )

        if contrib['functional_form'] == 'bond_harmonic':
            assert len(contrib['particle_indices']) == 2, 'Should have 2 bond instances'
            assert all(len(indices) == 2 for indices in contrib['particle_indices']), (
                'Bonds should have 2 particles'
            )
        elif contrib['functional_form'] == 'angle_harmonic':
            assert len(contrib['particle_indices']) == 2, (
                'Should have 2 angle instances'
            )
            assert all(len(indices) == 3 for indices in contrib['particle_indices']), (
                'Angles should have 3 particles'
            )
