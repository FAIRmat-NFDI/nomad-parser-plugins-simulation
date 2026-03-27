import numpy as np
from nomad_simulations.schema_packages.atoms_state import AtomsState, CGBeadState
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.utils.mdparserutils import (
    MDParser,
    particle_state_payloads_from_labels,
    particle_states_from_labels,
)


def test_particle_states_from_labels_promotes_atomic_labels():
    particle_states = particle_states_from_labels(['H', 'O'])

    assert [type(state) for state in particle_states] == [AtomsState, AtomsState]
    assert [state.chemical_symbol for state in particle_states] == ['H', 'O']
    assert [state.label for state in particle_states] == ['H', 'O']


def test_particle_states_from_labels_falls_back_to_cg():
    particle_states = particle_states_from_labels(['B1', 'monomer'])

    assert [type(state) for state in particle_states] == [CGBeadState, CGBeadState]
    assert [state.bead_symbol for state in particle_states] == ['B1', 'monomer']
    assert [state.label for state in particle_states] == ['B1', 'monomer']


def test_particle_state_payloads_from_labels_falls_back_to_cg_for_missing_labels():
    payloads = particle_state_payloads_from_labels(['H', None, ''])

    assert payloads == [
        {
            'm_def': CGBeadState.m_def.qualified_name(),
            'label': 'H',
            'bead_symbol': 'H',
        },
        {'m_def': CGBeadState.m_def.qualified_name(), 'label': None},
        {'m_def': CGBeadState.m_def.qualified_name(), 'label': ''},
    ]


def test_particle_state_payloads_from_numpy_array_labels():
    payloads = particle_state_payloads_from_labels(np.array(['H', 'O']))

    assert payloads == [
        {
            'm_def': AtomsState.m_def.qualified_name(),
            'chemical_symbol': 'H',
            'label': 'H',
        },
        {
            'm_def': AtomsState.m_def.qualified_name(),
            'chemical_symbol': 'O',
            'label': 'O',
        },
    ]


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
