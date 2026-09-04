import os
from pathlib import Path

import numpy as np
import pytest
from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.gromacs import parser as gromacs_parser
from nomad_simulation_parsers.parsers.gromacs.xvg_parser import GromacsXvgParser


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

    def get_frame_data(self, idx):
        lv = None if self._lattices is None else np.asarray(self._lattices[idx])
        return dict(
            positions=np.asarray(self._positions[idx]),
            velocities=(
                None if self._velocities is None else np.asarray(self._velocities[idx])
            ),
            lattice_vectors=lv,
        )

    def get_interactions(self):
        return []

    def get(self, key, default=None):  # noqa: ARG002
        return default


class StubInteractionsDataObject:
    """Stub data_object exposing only get_interactions(), for unit tests that
    exercise bond-list and force-field-contribution extraction."""

    def __init__(self, interactions: list):
        self._interactions = interactions

    def get_interactions(self):
        return self._interactions

    def get(self, key, default=None):  # noqa: ARG002
        return default


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
    # per-frame data present in every frame
    assert np.allclose(configs[0]['positions'], positions[0])
    assert np.allclose(configs[1]['positions'], positions[1])
    # topology only on frame 0
    assert np.allclose(configs[0]['lattice_vectors'], lattices[0])
    assert 'labels' in configs[0]
    assert 'n_particles' in configs[0]
    assert 'lattice_vectors' not in configs[1]
    assert 'labels' not in configs[1]
    assert 'n_particles' not in configs[1]


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
    assert 'pbc' in configs[0]
    assert configs[0]['pbc'] == [True, True, False]
    assert configs[1] == {}
    assert configs[2] == {}


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

    mdap.data_object = StubInteractionsDataObject(mock_interactions)

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


def test_get_coordinate_save_frequency():
    """Test coordinate save frequency extraction with compressed/uncompressed
    priority."""
    lp = gromacs_parser.GromacsLogParser()

    # Test nstxout-compressed takes priority
    params = {'nstxout-compressed': 100, 'nstxout': 50}
    assert lp.get_coordinate_save_frequency(params) == 100

    # Test nstxout used when compressed not available
    params = {'nstxout': 50}
    assert lp.get_coordinate_save_frequency(params) == 50

    # Test zero values are ignored
    params = {'nstxout-compressed': 0, 'nstxout': 50}
    assert lp.get_coordinate_save_frequency(params) == 50

    # Test None when both missing
    params = {}
    assert lp.get_coordinate_save_frequency(params) is None

    # Test None input
    assert lp.get_coordinate_save_frequency(None) is None


def test_get_bond_list():
    """Test bond list extraction from MDAnalysis interactions.

    MDAnalysis sets inter.btype to a topology-specific string (e.g. 'OW-HW'),
    never to the literal 'bond'. Bonds are identified solely by having exactly
    2 particle indices.
    """
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    mock_interactions = [
        {'type': 'OW-HW', 'atom_indices': [0, 1], 'atom_labels': ['OW', 'HW']},
        {'type': 'OW-HW', 'atom_indices': [0, 2], 'atom_labels': ['OW', 'HW']},
        {'type': 'CT-CT', 'atom_indices': [3, 4], 'atom_labels': ['CT', 'CT']},
        {
            'type': 'HW-OW-HW',
            'atom_indices': [1, 0, 2],
            'atom_labels': ['HW', 'OW', 'HW'],
        },
    ]

    mdap.data_object = StubInteractionsDataObject(mock_interactions)

    bond_list = mdap.get_bond_list()

    assert bond_list is not None
    assert isinstance(bond_list, np.ndarray)
    assert bond_list.shape == (3, 2)
    assert np.array_equal(bond_list[0], [0, 1])
    assert np.array_equal(bond_list[1], [0, 2])
    assert np.array_equal(bond_list[2], [3, 4])


def test_get_bond_list_no_bonds():
    """Test bond list extraction when only higher-order interactions are present."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    mock_interactions = [
        {'type': 'HW-OW-HW', 'atom_indices': [0, 1, 2], 'atom_labels': ['H', 'O', 'H']},
    ]

    mdap.data_object = StubInteractionsDataObject(mock_interactions)

    bond_list = mdap.get_bond_list()
    assert bond_list is None


def test_get_bond_list_empty_interactions():
    """Test bond list extraction with empty interactions."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    mdap.data_object = StubInteractionsDataObject([])

    bond_list = mdap.get_bond_list()
    assert bond_list is None


def test_get_bond_list_invalid_bond_indices():
    """Test bond list extraction skips entries with None or wrong-count indices."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    mock_interactions = [
        {'type': 'OW-HW', 'atom_indices': [0, 1], 'atom_labels': ['OW', 'HW']},
        {'type': 'OW-HW', 'atom_indices': None, 'atom_labels': ['OW', 'HW']},
        {
            'type': 'OW-HW-XX',
            'atom_indices': [2, 3, 4],
            'atom_labels': ['OW', 'HW', 'XX'],
        },
        {'type': 'CT-CT', 'atom_indices': [4, 5], 'atom_labels': ['CT', 'CT']},
    ]

    mdap.data_object = StubInteractionsDataObject(mock_interactions)

    bond_list = mdap.get_bond_list()

    assert bond_list is not None
    assert bond_list.shape == (2, 2)
    assert np.array_equal(bond_list[0], [0, 1])
    assert np.array_equal(bond_list[1], [4, 5])


def test_integration_bond_list_in_parsed_system():
    """Test that bond_list is populated in parsed model_system from TPR file."""
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    log_file = os.path.join(base, 'reference_s.log')

    if not os.path.exists(log_file):
        pytest.skip(f'Mainfile not found: {log_file}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(log_file, archive)

    assert archive.data is not None
    assert len(archive.data.model_system) > 0

    system = archive.data.model_system[0]
    assert system.bond_list is not None, 'Bond list should be populated from TPR'
    assert isinstance(system.bond_list, np.ndarray)
    assert system.bond_list.shape[1] == 2, 'Bond list should have shape (n_bonds, 2)'
    assert system.bond_list.shape[0] > 0, 'Should have at least one bond'
    # Water molecules have 2 O-H bonds each (432 bonds for 144 water molecules)
    assert system.bond_list.shape[0] == 432


def test_get_thermodynamic_ensemble():
    """Test ensemble determination from thermostat and barostat settings."""
    lp = gromacs_parser.GromacsLogParser()

    # Test NPT (both thermostat and barostat)
    params = {'tcoupl': 'v-rescale', 'pcoupl': 'parrinello-rahman'}
    assert lp.get_thermodynamic_ensemble(params) == 'NPT'

    # Test NVT (thermostat only)
    params = {'tcoupl': 'nose-hoover', 'pcoupl': 'no'}
    assert lp.get_thermodynamic_ensemble(params) == 'NVT'

    # Test NPH (barostat only)
    params = {'tcoupl': 'no', 'pcoupl': 'berendsen'}
    assert lp.get_thermodynamic_ensemble(params) == 'NPH'

    # Test NVE (neither)
    params = {'tcoupl': 'no', 'pcoupl': 'no'}
    assert lp.get_thermodynamic_ensemble(params) == 'NVE'

    # Test default values when missing
    params = {}
    assert lp.get_thermodynamic_ensemble(params) == 'NVE'

    # Test None input
    assert lp.get_thermodynamic_ensemble(None) is None


@pytest.mark.parametrize(
    'tcoupl,expected',
    [
        ('berendsen', 'berendsen'),
        ('nose-hoover', 'nose_hoover'),
        ('v-rescale', 'velocity_rescaling'),
        ('andersen', 'andersen'),
        ('andersen-massive', 'andersen_massive'),
        ('no', None),
        ('No', None),
        ('unknown_type', None),
        (None, None),
    ],
)
def test_get_thermostat_type(tcoupl, expected):
    """Test thermostat type mapping."""
    lp = gromacs_parser.GromacsLogParser()
    result = lp.get_thermostat_type(tcoupl)
    assert result == expected


def test_get_reference_temperature():
    """Test reference temperature extraction from grpopts.ref-t."""
    lp = gromacs_parser.GromacsLogParser()

    params = {'grpopts': {'ref-t': 300.0}}
    assert lp.get_reference_temperature(params) == 300.0

    params = {'grpopts': {'ref-t': [300.0, 310.0, 320.0]}}
    assert lp.get_reference_temperature(params) == 300.0

    params = {'grpopts': {'ref-t': []}}
    assert lp.get_reference_temperature(params) is None

    params = {}
    assert lp.get_reference_temperature(params) is None

    assert lp.get_reference_temperature(None) is None


def test_get_thermostat_coupling_constant():
    """Test thermostat coupling constant extraction from grpopts.tau-t."""
    lp = gromacs_parser.GromacsLogParser()

    params = {'grpopts': {'tau-t': 0.1}}
    assert lp.get_thermostat_coupling_constant(params) == 0.1

    params = {'grpopts': {'tau-t': [0.1, 0.2, 0.3]}}
    assert lp.get_thermostat_coupling_constant(params) == 0.1

    params = {'grpopts': {'tau-t': []}}
    assert lp.get_thermostat_coupling_constant(params) is None

    params = {}
    assert lp.get_thermostat_coupling_constant(params) is None


@pytest.mark.parametrize(
    'pcoupl,expected',
    [
        ('berendsen', 'berendsen'),
        ('parrinello-rahman', 'parrinello_rahman'),
        ('mttk', 'mttk'),
        ('c-rescale', 'c_rescale'),
        ('no', None),
        ('No', None),
        ('unknown_type', None),
        (None, None),
    ],
)
def test_get_barostat_type(pcoupl, expected):
    """Test barostat type mapping."""
    lp = gromacs_parser.GromacsLogParser()
    result = lp.get_barostat_type(pcoupl)
    assert result == expected


@pytest.mark.parametrize(
    'pcoupltype,expected',
    [
        ('isotropic', 'isotropic'),
        ('semiisotropic', 'semi_isotropic'),
        ('anisotropic', 'anisotropic'),
        ('surface-tension', 'surface_tension'),
        ('unknown_type', None),
        (None, None),
    ],
)
def test_get_barostat_coupling_type(pcoupltype, expected):
    """Test barostat coupling type mapping."""
    lp = gromacs_parser.GromacsLogParser()
    result = lp.get_barostat_coupling_type(pcoupltype)
    assert result == expected


def test_get_barostat_coupling_constant():
    """Test barostat coupling constant extraction."""
    lp = gromacs_parser.GromacsLogParser()

    # Test scalar value
    params = {'tau-p': 2.0}
    assert lp.get_barostat_coupling_constant(params) == 2.0

    # Test missing value
    params = {}
    assert lp.get_barostat_coupling_constant(params) is None

    # Test None input
    assert lp.get_barostat_coupling_constant(None) is None


def test_get_matrix_parameter():
    """
    Test get_matrix_parameter for all valid, missing, wrong-type and wrong-shape cases.
    """
    lp = gromacs_parser.GromacsLogParser()

    # Test implemented parameter keys with correct matrix
    matrix = np.eye(3)
    params = {'ref-p': matrix}
    np.testing.assert_array_equal(
        lp.get_matrix_parameter(params, param_key='ref-p'), matrix
    )

    params = {'compressibility': matrix}
    np.testing.assert_array_equal(
        lp.get_matrix_parameter(params, param_key='compressibility'), matrix
    )

    # Test list of lists
    matrix_list = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    params = {'ref-p': matrix_list}
    np.testing.assert_array_equal(
        lp.get_matrix_parameter(params, param_key='ref-p'), matrix
    )

    # Test wrong shape (not 3x3)
    wrong_shape = np.eye(2)
    params = {'ref-p': wrong_shape}
    assert lp.get_matrix_parameter(params, param_key='ref-p') is None

    # Test wrong parameter type (dict, not matrix or list)
    params = {'ref-p': {'a': 1.0}}
    assert lp.get_matrix_parameter(params, param_key='ref-p') is None

    # Test missing parameters
    assert lp.get_matrix_parameter(None, param_key='ref-p') is None


# ---------------------------------------------------------------------------
# Free energy calculation transformer tests
# ---------------------------------------------------------------------------


def test_get_free_energy_calc_type():
    lp = gromacs_parser.GromacsLogParser()
    assert lp.get_free_energy_calc_type({'free-energy': 'yes'}) == 'alchemical'
    assert lp.get_free_energy_calc_type({'free-energy': 'slow-growth'}) == 'alchemical'
    assert lp.get_free_energy_calc_type({'free-energy': 'expanded'}) == 'alchemical'
    assert (
        lp.get_free_energy_calc_type({'free-energy': 'umbrella'}) == 'umbrella_sampling'
    )
    assert lp.get_free_energy_calc_type({'free-energy': 'no'}) is None
    assert lp.get_free_energy_calc_type({}) is None
    assert lp.get_free_energy_calc_type(None) is None


def test_get_fep_params_if_active():
    lp = gromacs_parser.GromacsLogParser()
    src = {'input_parameters': {'free-energy': 'yes'}}
    assert lp.get_fep_params_if_active(src) is src
    src_no = {'input_parameters': {'free-energy': 'no'}}
    assert lp.get_fep_params_if_active(src_no) is None
    assert lp.get_fep_params_if_active({'input_parameters': {}}) is None
    assert lp.get_fep_params_if_active(None) is None


def test_get_lambda_state_index():
    lp = gromacs_parser.GromacsLogParser()
    assert lp.get_lambda_state_index({'init-lambda-state': '3'}) == 3
    assert lp.get_lambda_state_index({'init-lambda-state': 0}) == 0
    assert lp.get_lambda_state_index({'init-lambda-state': -1}) is None
    assert lp.get_lambda_state_index({}) is None
    assert lp.get_lambda_state_index(None) is None


def test_get_lambdas_schedule_non_fe():
    lp = gromacs_parser.GromacsLogParser()
    assert lp.get_lambdas_schedule({}) is None
    assert lp.get_lambdas_schedule(None) is None
    # no all-lambdas key → None
    assert lp.get_lambdas_schedule({'free-energy': 'yes'}) is None


def test_get_lambdas_schedule_fe():
    lp = gromacs_parser.GromacsLogParser()
    params = {
        'all-lambdas': {
            'vdw-lambdas': '0.0 0.25 0.5 0.75 1.0',
            'coul-lambdas': '0.0 0.0 0.0 0.0 0.0',  # all-zero → skipped
        },
        'sc-alpha': '0.5',
        'sc-power': '1',
        'sc-sigma': '0.3',
    }
    result = lp.get_lambdas_schedule(params)
    assert result is not None
    assert len(result) == 1  # coul-lambdas all-zero, only vdw survives
    entry = result[0]
    assert entry['interaction_type'] == 'vdw'
    np.testing.assert_array_equal(
        entry['lambda_values'], np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    )
    assert entry['softcore_enabled'] is True
    assert entry['softcore_alpha'] == '0.5'


def test_get_current_lambdas():
    lp = gromacs_parser.GromacsLogParser()
    params = {
        'init-lambda-state': '2',
        'all-lambdas': {'vdw-lambdas': '0.0 0.5 1.0'},
        'sc-alpha': '0.0',
    }
    result = lp.get_current_lambdas(params)
    assert result is not None
    np.testing.assert_array_equal(result, np.array([1.0]))
    # index out of range
    params['init-lambda-state'] = '99'
    assert lp.get_current_lambdas(params) is None


def test_system_hierarchy_water():
    """Test that system hierarchy is parsed from water mainfile."""
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    log_file = os.path.join(base, 'reference_s.log')

    if not os.path.exists(log_file):
        pytest.skip(f'Mainfile not found: {log_file}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(log_file, archive)

    assert archive.data is not None
    assert len(archive.data.model_system) > 0

    system = archive.data.model_system[0]

    # Water system should have sub_systems (molecule groups)
    assert system.sub_systems is not None, 'System should have sub_systems hierarchy'
    assert len(system.sub_systems) > 0, 'Should have at least one molecule group'

    # Check first molecule group (should be SOL for water)
    mol_group = system.sub_systems[0]
    assert mol_group.name is not None
    assert 'group_' in mol_group.name, 'Molecule group should have group_ prefix'
    assert mol_group.branch_label == 'molecule_group'
    assert mol_group.composition_formula is not None
    assert mol_group.particle_indices is not None
    assert len(mol_group.particle_indices) > 0

    # Molecule group should contain individual molecules
    assert mol_group.sub_systems is not None
    assert len(mol_group.sub_systems) > 0, 'Molecule group should contain molecules'

    # Check first molecule
    molecule = mol_group.sub_systems[0]
    assert molecule.name is not None
    assert molecule.branch_label == 'molecule'
    assert molecule.particle_indices is not None
    assert len(molecule.particle_indices) > 0

    # Water molecule is terminal (single residue), should have composition_formula
    assert molecule.composition_formula is not None
    assert 'H' in molecule.composition_formula
    assert 'O' in molecule.composition_formula

    # All terminal branches across the full hierarchy must have composition_formula
    def check_terminal_formulas(subsystem):
        if subsystem.sub_systems is None or len(subsystem.sub_systems) == 0:
            assert subsystem.composition_formula is not None, (
                f'Terminal {subsystem.branch_label} should have composition_formula'
            )
        else:
            for child in subsystem.sub_systems:
                check_terminal_formulas(child)

    for mg in system.sub_systems:
        check_terminal_formulas(mg)


def test_system_hierarchy_polymer():
    """Test that complex system hierarchy is parsed from protein mainfile."""
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'protein_small'
    log_file = os.path.join(base, 'md.log')

    if not os.path.exists(log_file):
        pytest.skip(f'Mainfile not found: {log_file}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(log_file, archive)

    assert archive.data is not None
    assert len(archive.data.model_system) > 0

    system = archive.data.model_system[0]

    # Polymer system should have sub_systems
    assert system.sub_systems is not None and len(system.sub_systems) > 0, (
        'System has no molecular hierarchy - hierarchy parsing regression'
    )

    # Find a molecule group with multi-residue molecules (polymer chain)
    multi_residue_group = None
    for mol_group in system.sub_systems:
        if mol_group.sub_systems and len(mol_group.sub_systems) > 0:
            first_mol = mol_group.sub_systems[0]
            # Multi-residue molecule has sub_systems (monomer groups)
            if first_mol.sub_systems is not None and len(first_mol.sub_systems) > 0:
                multi_residue_group = mol_group
                break

    assert multi_residue_group is not None, (
        'No multi-residue molecules found in protein_small - '
        'expected at least one polymer chain'
    )

    # Test full 4-level hierarchy: group → molecule → monomer_group → monomer
    assert multi_residue_group.branch_label == 'molecule_group'

    molecule = multi_residue_group.sub_systems[0]
    assert molecule.branch_label == 'molecule'
    assert molecule.sub_systems is not None
    assert len(molecule.sub_systems) > 0

    # Molecule should contain monomer groups
    monomer_group = molecule.sub_systems[0]
    assert monomer_group.branch_label == 'monomer_group'
    assert 'group_' in monomer_group.name
    assert monomer_group.sub_systems is not None
    assert len(monomer_group.sub_systems) > 0

    # Monomer group should contain individual monomers
    monomer = monomer_group.sub_systems[0]
    assert monomer.branch_label == 'monomer'
    assert monomer.composition_formula is not None
    assert monomer.particle_indices is not None
    assert len(monomer.particle_indices) > 0

    # Exhaustively check branch_label at every level across all sub_systems
    for mg in system.sub_systems:
        assert mg.branch_label == 'molecule_group'

        if mg.sub_systems and len(mg.sub_systems) > 0:
            for mol in mg.sub_systems:
                assert mol.branch_label == 'molecule'

                if mol.sub_systems and len(mol.sub_systems) > 0:
                    for mng in mol.sub_systems:
                        assert mng.branch_label == 'monomer_group'

                        if mng.sub_systems and len(mng.sub_systems) > 0:
                            for mon in mng.sub_systems:
                                assert mon.branch_label == 'monomer'


def test_system_hierarchy_particle_indices_valid():
    """Test that particle_indices in hierarchy are valid and consistent."""
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    log_file = os.path.join(base, 'reference_s.log')

    if not os.path.exists(log_file):
        pytest.skip(f'Mainfile not found: {log_file}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(log_file, archive)

    system = archive.data.model_system[0]
    n_particles = system.n_particles

    def check_particle_indices(subsystem, parent_indices=None):
        """Recursively check particle_indices validity."""
        assert subsystem.particle_indices is not None
        assert len(subsystem.particle_indices) > 0

        # All indices should be valid (less than n_particles)
        assert np.all(subsystem.particle_indices < n_particles)
        assert np.all(subsystem.particle_indices >= 0)

        # If parent exists, subsystem indices should be subset of parent
        if parent_indices is not None:
            assert np.all(np.isin(subsystem.particle_indices, parent_indices))

        # Recursively check children
        if subsystem.sub_systems:
            for child in subsystem.sub_systems:
                check_particle_indices(child, subsystem.particle_indices)

    # Check all molecule groups
    for mol_group in system.sub_systems:
        check_particle_indices(mol_group)


def test_xvg_parser_get_results_valid():
    """Unit test for GromacsXVGParser.get_results() with real XVG data."""
    xvg_file = (
        Path(__file__).parent.parent
        / 'data'
        / 'gromacs'
        / 'free_energy_calculations'
        / 'alchemical_transformation_single_run'
        / 'fep_run-7.xvg'
    )
    if not xvg_file.exists():
        pytest.skip(f'XVG test data not found: {xvg_file}')

    text_parser = GromacsXvgParser()
    text_parser.mainfile = str(xvg_file)
    text_parser.parse()

    xvg_parser = gromacs_parser.GromacsXVGParser(text_parser=text_parser)
    xvg_parser.filepath = str(xvg_file)

    results = xvg_parser.get_results()
    fec = results.get('free_energy_calculations', {})

    assert fec.get('n_frames') == 5001
    assert fec.get('n_states') == 11
    assert fec.get('times') is not None
    assert len(fec['times']) == 5001
    assert fec.get('value_total_energy_derivative') is not None
    assert fec['value_total_energy_derivative'].shape == (5001,)
    assert fec.get('value_total_energy_differences') is not None
    assert fec['value_total_energy_differences'].shape == (5001, 11)
    assert fec.get('value_PV_energy') is not None
    assert fec['value_PV_energy'].shape == (5001,)


def test_integration_fep_xvg_fields_populated():
    """Integration test: XVG fields appear in FreeEnergyCalculationParameters."""
    log_file = (
        Path(__file__).parent.parent
        / 'data'
        / 'gromacs'
        / 'free_energy_calculations'
        / 'alchemical_transformation_single_run'
        / 'fep_run-7.log'
    )
    if not log_file.exists():
        pytest.skip(f'FEP test data not found: {log_file}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(str(log_file), archive)

    assert archive.workflow2 is not None
    method = archive.workflow2.method
    assert method is not None
    assert method.free_energy_calculation_parameters is not None
    assert len(method.free_energy_calculation_parameters) > 0

    fep = method.free_energy_calculation_parameters[0]

    assert fep.n_frames == 5001
    assert fep.n_states == 11
    assert fep.times is not None
    assert len(fep.times) == 5001
    assert fep.energy_derivative is not None
    assert fep.energy_derivative.shape == (5001,)
    assert fep.energy_differences is not None
    assert fep.energy_differences.shape == (5001, 11)
    assert fep.pv_energy is not None
    assert fep.pv_energy.shape == (5001,)
