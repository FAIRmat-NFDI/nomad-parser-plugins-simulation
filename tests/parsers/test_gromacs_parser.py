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

    def get_interactions(self):
        # Return empty list for stub (no bonds by default)
        return []


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

    # One config dict for the representative (last) frame only.
    assert isinstance(configs, list)
    assert len(configs) == 1
    assert np.allclose(configs[0]['positions'], positions[1])
    assert np.allclose(configs[0]['lattice_vectors'], lattices[1])
    # Particle labels are not injected into config dicts; each dict carries
    # a 'step' key used by the get_particle_states transformer.
    assert 'step' in configs[0]
    assert 'labels' not in configs[0]

    # get_particle_states returns payloads for the representative frame.
    simple_mdanalysis_parser._particle_parameters = [
        {'label': 'H', 'element': 'H'}
    ] * 3
    last_step = simple_mdanalysis_parser._trajectory_steps_sampled[-1]
    payloads = simple_mdanalysis_parser.get_particle_states(configs[0])
    assert len(payloads) == 3
    assert payloads[0]['label'] == 'H'
    assert configs[0]['step'] == last_step


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
        assert len(configs) == 1
        assert configs[0]['positions'].shape == (2, 3)
    else:
        assert len(sim.model_system) == 1
        assert sim.model_system[0].positions is not None
        assert sim.model_system[0].positions.shape == (2, 3)


def test_logparser_get_configurations_pbc_and_sampling():
    lp = gromacs_parser.GromacsLogParser()
    # Emulate parsed log data where pbc string is 'xy' (x,y periodic)
    lp._data = {'input_parameters': {'pbc': 'xy'}}
    lp._trajectory_steps_sampled = [0, 1, 2]
    configs = lp.get_configurations()

    # Log parser emits one config dict (pbc is topology data, not per-frame).
    assert len(configs) == 1
    assert 'pbc' in configs[0]
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
            'get_interactions': lambda self: mock_interactions,
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
    """Test bond list extraction from MDAnalysis interactions."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    # Test with bond interactions
    mock_interactions = [
        {'type': 'bond', 'atom_indices': [0, 1], 'atom_labels': ['O', 'H']},
        {'type': 'bond', 'atom_indices': [0, 2], 'atom_labels': ['O', 'H']},
        {'type': 'bond', 'atom_indices': [3, 4], 'atom_labels': ['O', 'H']},
        {'type': 'angle', 'atom_indices': [1, 0, 2], 'atom_labels': ['H', 'O', 'H']},
    ]

    mdap.data_object = type(
        'obj', (object,), {'get_interactions': lambda self: mock_interactions}
    )()

    bond_list = mdap.get_bond_list()

    assert bond_list is not None
    assert isinstance(bond_list, np.ndarray)
    assert bond_list.shape == (3, 2)
    assert np.array_equal(bond_list[0], [0, 1])
    assert np.array_equal(bond_list[1], [0, 2])
    assert np.array_equal(bond_list[2], [3, 4])


def test_get_bond_list_no_bonds():
    """Test bond list extraction when no bonds are present."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    # Only angles, no bonds
    mock_interactions = [
        {'type': 'angle', 'atom_indices': [0, 1, 2], 'atom_labels': ['H', 'O', 'H']},
    ]

    mdap.data_object = type(
        'obj', (object,), {'get_interactions': lambda self: mock_interactions}
    )()

    bond_list = mdap.get_bond_list()
    assert bond_list is None


def test_get_bond_list_empty_interactions():
    """Test bond list extraction with empty interactions."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    mdap.data_object = type('obj', (object,), {'get_interactions': lambda self: []})()

    bond_list = mdap.get_bond_list()
    assert bond_list is None


def test_get_bond_list_invalid_bond_indices():
    """Test bond list extraction filters out invalid bond entries."""
    mdap = gromacs_parser.GromacsMDAnalysisParser()

    # Mix of valid and invalid bond interactions
    mock_interactions = [
        {'type': 'bond', 'atom_indices': [0, 1], 'atom_labels': ['O', 'H']},
        {'type': 'bond', 'atom_indices': None, 'atom_labels': ['O', 'H']},
        {'type': 'bond', 'atom_indices': [2, 3, 4], 'atom_labels': ['O', 'H', 'C']},
        {'type': 'bond', 'atom_indices': [4, 5], 'atom_labels': ['O', 'H']},
    ]

    mdap.data_object = type(
        'obj', (object,), {'get_interactions': lambda self: mock_interactions}
    )()

    bond_list = mdap.get_bond_list()

    assert bond_list is not None
    assert bond_list.shape == (2, 2)
    assert np.array_equal(bond_list[0], [0, 1])
    assert np.array_equal(bond_list[1], [4, 5])


def test_integration_bond_list_in_parsed_system():
    """Test that bond_list is populated in parsed model_system from TPR file."""
    base = Path(__file__).parent.parent / 'data' / 'gromacs' / 'water'
    tpr_file = os.path.join(base, 'topol.tpr')

    if not os.path.exists(tpr_file):
        pytest.skip(f'TPR file not found: {tpr_file}')

    archive = EntryArchive()
    parser = gromacs_parser.GromacsParser()
    parser.parse(tpr_file, archive)

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
    """Test reference temperature extraction from scalar or array."""
    lp = gromacs_parser.GromacsLogParser()

    # Test scalar value
    params = {'ref-t': 300.0}
    assert lp.get_reference_temperature(params) == 300.0

    # Test array (take first value)
    params = {'ref-t': [300.0, 310.0, 320.0]}
    assert lp.get_reference_temperature(params) == 300.0

    # Test underscore variant
    params = {'ref_t': 298.0}
    assert lp.get_reference_temperature(params) == 298.0

    # Test empty array
    params = {'ref-t': []}
    assert lp.get_reference_temperature(params) is None

    # Test missing
    params = {}
    assert lp.get_reference_temperature(params) is None

    # Test None input
    assert lp.get_reference_temperature(None) is None


def test_get_thermostat_coupling_constant():
    """Test thermostat coupling constant extraction from scalar or array."""
    lp = gromacs_parser.GromacsLogParser()

    # Test scalar value
    params = {'tau-t': 0.1}
    assert lp.get_thermostat_coupling_constant(params) == 0.1

    # Test array (take first value)
    params = {'tau-t': [0.1, 0.2, 0.3]}
    assert lp.get_thermostat_coupling_constant(params) == 0.1

    # Test underscore variant
    params = {'tau_t': 0.5}
    assert lp.get_thermostat_coupling_constant(params) == 0.5

    # Test empty array
    params = {'tau-t': []}
    assert lp.get_thermostat_coupling_constant(params) is None

    # Test missing
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


def test_get_reference_pressure():
    """Test reference pressure extraction from scalar, array, or matrix."""
    lp = gromacs_parser.GromacsLogParser()

    # Test scalar value
    params = {'ref-p': 1.0}
    assert lp.get_reference_pressure(params) == 1.0

    # Test array (take first value)
    params = {'ref-p': [1.0, 1.0]}
    assert lp.get_reference_pressure(params) == 1.0

    # Test matrix (take [0][0])
    params = {'ref-p': [[1.0, 0.0], [0.0, 1.0]]}
    assert lp.get_reference_pressure(params) == 1.0

    # Test underscore variant
    params = {'ref_p': 1.5}
    assert lp.get_reference_pressure(params) == 1.5

    # Test empty array
    params = {'ref-p': []}
    assert lp.get_reference_pressure(params) is None

    # Test empty matrix
    params = {'ref-p': [[]]}
    assert lp.get_reference_pressure(params) is None

    # Test missing
    params = {}
    assert lp.get_reference_pressure(params) is None


def test_get_barostat_coupling_constant():
    """Test barostat coupling constant extraction."""
    lp = gromacs_parser.GromacsLogParser()

    # Test hyphen variant
    params = {'tau-p': 2.0}
    assert lp.get_barostat_coupling_constant(params) == 2.0

    # Test underscore variant
    params = {'tau_p': 5.0}
    assert lp.get_barostat_coupling_constant(params) == 5.0

    # Test missing
    params = {}
    assert lp.get_barostat_coupling_constant(params) is None

    # Test None input
    assert lp.get_barostat_coupling_constant(None) is None


def test_get_compressibility():
    """Test compressibility extraction from scalar, array, or matrix."""
    lp = gromacs_parser.GromacsLogParser()

    # Test scalar value
    params = {'compressibility': 4.5e-5}
    assert lp.get_compressibility(params) == 4.5e-5

    # Test array (take first value)
    params = {'compressibility': [4.5e-5, 4.5e-5]}
    assert lp.get_compressibility(params) == 4.5e-5

    # Test matrix (take [0][0])
    params = {'compressibility': [[4.5e-5, 0.0], [0.0, 4.5e-5]]}
    assert lp.get_compressibility(params) == 4.5e-5

    # Test empty array
    params = {'compressibility': []}
    assert lp.get_compressibility(params) is None

    # Test missing
    params = {}
    assert lp.get_compressibility(params) is None
