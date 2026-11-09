import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.gaussian.parser import OutParser
from tests.parsers.common import approx, assert_approx


@pytest.mark.unit
class TestGaussianMapping:
    def test_maps_program_metadata(self):
        program = OutParser().get_program_data(
            {'program': ['09', 'EM64L-G09RevD.01', '24-Apr-2013']}
        )

        assert program == {
            'name': 'Gaussian',
            'version': '09 EM64L-G09RevD.01 24-Apr-2013',
        }

    def test_maps_only_dft_routes(self):
        methods = OutParser().get_dft_method(
            {
                'run': [
                    {'x_gaussian_settings_corrected': '#B3LYP/aug-cc-pVTZ'},
                    {'x_gaussian_settings_corrected': '#CCSD(T)/6-31G*'},
                ]
            }
        )

        assert methods == [{'name': 'DFT', 'xc': {'functional_key': 'B3LYP'}}]

    def test_maps_atomic_numbers_and_positions(self):
        systems = OutParser().get_atoms(
            {
                'run': [
                    {
                        'system': [
                            {
                                'standard_orientation': np.array(
                                    [
                                        [1, 13, 0, 0.0, 0.0, 1.25],
                                        [2, 8, 0, 0.0, 0.0, -1.25],
                                    ]
                                )
                            }
                        ]
                    }
                ]
            }
        )

        assert [state['atomic_number'] for state in systems[0]['particle_states']] == [
            13,
            8,
        ]
        assert [
            state['chemical_symbol'] for state in systems[0]['particle_states']
        ] == ['Al', 'O']
        assert_approx(
            systems[0]['positions'].to('angstrom').magnitude,
            [[0.0, 0.0, 1.25], [0.0, 0.0, -1.25]],
        )

    def test_maps_particle_identity_only_to_topology_frame(self):
        orientation = np.array([[1, 13, 0, 0.0, 0.0, 0.0]])
        systems = OutParser().get_atoms(
            {
                'run': [
                    {
                        'system': [
                            {'standard_orientation': orientation},
                            {'standard_orientation': orientation.copy()},
                        ]
                    }
                ]
            }
        )

        assert systems[0]['particle_states'] == [
            {'atomic_number': 13, 'chemical_symbol': 'Al'}
        ]
        assert 'particle_states' not in systems[1]

    def test_maps_energy_and_forces(self):
        outputs = OutParser().get_outputs(
            {
                'run': [
                    {
                        'system': [
                            {
                                'scf_iteration_final': {
                                    'energy_total': -5.0 * ureg.hartree
                                },
                                'forces': np.array([[1.0, 2.0, 3.0]]),
                            }
                        ]
                    }
                ]
            }
        )

        assert outputs[0]['total_energies'][0]['value'].to(
            'hartree'
        ).magnitude == approx(-5.0)
        assert_approx(outputs[0]['total_forces'][0]['value'], [[1.0, 2.0, 3.0]])
