import numpy as np
from nomad.units import ureg
from nomad_simulations.schema_packages.atoms_state import AtomsState, CGBeadState
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser


def test_particle_states_from_labels_promotes_atomic_labels():
    parameters = MDParser.particle_parameters_from_labels(['H', 'O'])
    particle_states = MDParser.particle_states_from_parameters(parameters)

    assert [type(state) for state in particle_states] == [AtomsState, AtomsState]
    assert [state.chemical_symbol for state in particle_states] == ['H', 'O']
    assert [state.label for state in particle_states] == ['H', 'O']


def test_particle_states_from_labels_falls_back_to_cg():
    parameters = MDParser.particle_parameters_from_labels(['B1', 'monomer'])
    particle_states = MDParser.particle_states_from_parameters(parameters)

    assert [type(state) for state in particle_states] == [CGBeadState, CGBeadState]
    assert [state.bead_symbol for state in particle_states] == ['B1', 'monomer']
    assert [state.label for state in particle_states] == ['B1', 'monomer']


def test_particle_state_payloads_from_labels_falls_back_to_cg_for_missing_labels():
    parameters = MDParser.particle_parameters_from_labels(['H', None, ''])
    payloads = MDParser.particle_state_payloads_from_parameters(parameters)

    assert payloads == [
        {
            'label': 'H',
            'bead_symbol': 'H',
        },
        {'label': None},
        {'label': ''},
    ]


def test_particle_state_payloads_from_numpy_array_labels():
    parameters = MDParser.particle_parameters_from_labels(np.array(['H', 'O']))
    payloads = MDParser.particle_state_payloads_from_parameters(parameters)

    assert payloads == [
        {
            'chemical_symbol': 'H',
            'label': 'H',
        },
        {
            'chemical_symbol': 'O',
            'label': 'O',
        },
    ]


def test_particle_states_from_parameters_supports_mixed_identity_and_mass():
    particle_states = MDParser.particle_states_from_parameters(
        [
            {'label': 'OW', 'chemical_symbol': 'O', 'mass': 15.999},
            {'label': 'WN1', 'mass': 50.0},
        ]
    )

    assert [type(state) for state in particle_states] == [AtomsState, CGBeadState]
    assert [state.get_label() for state in particle_states] == ['O', 'WN1']
    assert [state.label for state in particle_states] == ['OW', 'WN1']
    assert [state.mass.to(ureg.amu).magnitude for state in particle_states] == [
        15.999,
        50.0,
    ]


def test_element_like_label_remains_cg_without_explicit_chemical_symbol():
    [particle_state] = MDParser.particle_states_from_parameters(
        [{'label': 'C', 'mass': 50.0, 'charge': -0.5}]
    )

    assert type(particle_state) is CGBeadState
    assert particle_state.bead_symbol == 'C'
    assert particle_state.mass.to(ureg.amu).magnitude == 50.0
    assert particle_state.charge is None


def test_particle_states_ignore_invalid_mass_and_partial_charge():
    atom, bead = MDParser.particle_states_from_parameters(
        [
            {'label': 'OW', 'chemical_symbol': 'O', 'mass': np.nan, 'charge': -0.8},
            {'label': 'WN1', 'chemical_symbol': 'WN1', 'mass': 0, 'charge': 0.5},
        ]
    )

    assert atom.mass is None
    assert atom.charge == 0
    assert bead.mass is None
    assert bead.charge is None


def test_parse_trajectory_step_uses_atomic_particle_states():
    parser = MDParser()
    parser.trajectory_steps = [0]
    simulation = Simulation()

    parser.parse_trajectory_step(
        data={
            'step': 0,
            'positions': np.zeros((2, 3)),
            'lattice_vectors': np.eye(3),
            'periodic_boundary_conditions': [True, True, True],
            'labels': ['H', 'O'],
            'dimensions': 3,
        },
        simulation=simulation,
    )

    particle_states = simulation.model_system[0].particle_states
    assert [type(state) for state in particle_states] == [AtomsState, AtomsState]
    assert [state.chemical_symbol for state in particle_states] == ['H', 'O']


def test_parse_trajectory_step_uses_cg_particle_states():
    parser = MDParser()
    parser.trajectory_steps = [0]
    simulation = Simulation()

    parser.parse_trajectory_step(
        data={
            'step': 0,
            'positions': np.zeros((2, 3)),
            'lattice_vectors': np.eye(3),
            'periodic_boundary_conditions': [True, True, True],
            'labels': ['B1', 'monomer'],
            'dimensions': 3,
        },
        simulation=simulation,
    )

    particle_states = simulation.model_system[0].particle_states
    assert [type(state) for state in particle_states] == [CGBeadState, CGBeadState]
    assert [state.bead_symbol for state in particle_states] == ['B1', 'monomer']
