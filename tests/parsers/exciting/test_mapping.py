from unittest.mock import PropertyMock, patch

import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.outputs import Outputs, SCFSteps
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)

from nomad_simulation_parsers.parsers.exciting.parser import (
    BandstructureXMLParser,
    DosXMLParser,
    EigvalParser,
    InfoParser,
)


@pytest.mark.unit
class TestInfoParserMapping:
    def test_maps_xc_functionals_used_by_example_assertions(self):
        parser = InfoParser()
        assert parser.get_xc_functional_components(408) == [
            {'libxc': 'HYB_GGA_XC_HSE06'}
        ]
        assert parser.get_xc_functional_components(20) == [
            {'libxc': 'GGA_C_PBE'},
            {'libxc': 'GGA_X_PBE'},
        ]

    def test_maps_scf_convergence_quantities_with_units(self):
        parser = InfoParser()
        source = {
            'groundstate': {
                'scf_iteration': [
                    {
                        'energy_total': -2.0 * ureg.hartree,
                        'time_physical': 1.0 * ureg.second,
                        'x_exciting_energy_convergence': np.array([1e-4, 2e-4])
                        * ureg.hartree,
                        'x_exciting_effective_potential_convergence': np.array(
                            [3e-5, 4e-5]
                        )
                        * ureg.hartree,
                    },
                    {
                        'energy_total': -2.1 * ureg.hartree,
                        'time_physical': 3.0 * ureg.second,
                        'x_exciting_energy_convergence': np.array([5e-6, 6e-6])
                        * ureg.hartree,
                    },
                ]
            }
        }
        result = parser.get_scf_steps(source)
        np.testing.assert_allclose(result['durations'], [1.0, 2.0])
        assert len(result['energies_total']) == 2
        assert result['delta_energies_total'][-1].to(
            'hartree'
        ).magnitude == pytest.approx(5e-6)

    @pytest.mark.parametrize(
        ('scf_history', 'thresholds', 'expected_targets', 'expected_converged'),
        [
            pytest.param(
                # Total energy, RMS potential change, absolute charge change in a.u.
                [
                    (-470.41895633, 0, 0),
                    (-469.77639935, 16.3619, 0.00426985),
                    (-469.42507007, 12.0457, 0.00241697),
                    (-470.20883560, 0.086184, 0.00457412),
                    (-467.88686660, 0.107, 0.014833),
                    (-469.12499436, 0.00445151, 0.00796641),
                    (-469.21422156, 0.00310446, 0.00055505),
                    (-469.15655455, 0.000436029, 0.000364686),
                    (-469.16327297, 6.80138e-6, 4.16918e-5),
                    (-469.16317864, 3.39085e-6, 3.36169e-7),
                    (-469.16315587, 3.42312e-7, 2.16994e-7),
                    (-469.16316142, 7.43429e-8, 4.63388e-8),
                    (-469.16316077, 1.25638e-8, 5.26988e-9),
                    (-469.16316066, 3.71179e-9, 8.62576e-10),
                ],
                # Threshold values and residual column indices.
                {
                    'x_exciting_energy_convergence': (1e-6 * ureg.hartree, -1),
                    'x_exciting_effective_potential_convergence': (
                        1e-6 * ureg.hartree,
                        1,
                    ),
                    'x_exciting_charge_convergence': (1e-5 * ureg.elementary_charge, 2),
                },
                [
                    ('PotentialConvergenceTarget', True),
                    ('EnergyConvergenceTarget', True),
                    ('DensityConvergenceTarget', True),
                ],
                True,
                id='all-targets-met',
            ),
        ],
    )
    def test_scf_convergence(
        self, scf_history, thresholds, expected_targets, expected_converged
    ):
        # Add a column for the energy change between iterations
        scf_history = np.pad(
            scf_history,
            ((0, 0), (0, 1)),
            mode='constant',
            constant_values=0.0,
        )
        scf_history[1:, -1] = np.abs(np.diff(scf_history[:, 0], axis=0))

        iterations = []
        for step, iter_data in enumerate(scf_history):
            iteration = {
                'energy_total': iter_data[0] * ureg.hartree,
                'time_physical': (step + 1) * ureg.second,
            }
            if step > 0:
                for key, (threshold, icol) in thresholds.items():
                    iteration[key] = (
                        np.array([iter_data[icol], threshold.magnitude])
                        * threshold.units
                    )
            iterations.append(iteration)

        source = {'groundstate': {'scf_iteration': iterations}}
        parser = InfoParser()
        result = parser.get_scf_steps(source)

        workflow = SinglePoint(
            method=SinglePointMethod(
                convergence_targets=parser.get_single_point_convergence(source)
            )
        )
        archive = EntryArchive(
            data=Simulation(outputs=[Outputs(scf_steps=SCFSteps(**result))]),
            workflow2=workflow,
        )
        workflow.normalize(archive, get_logger(__name__))

        assert [
            (entry.convergence_target_ref.m_def.name, entry.is_reached)
            for entry in workflow.results.convergence
        ] == expected_targets
        assert workflow.results.is_converged is expected_converged

    def test_maps_smearing_kind(self):
        parser = InfoParser()
        assert (
            parser.get_smearing_kind('Extended linear tetrahedron method')
            == 'Tetrahedra-extended'
        )

    def test_maps_k_mesh(self):
        parser = InfoParser()
        source = {
            'initialization': {
                'kpoint_grid': [4, 4, 4],
                'x_exciting_number_kpoints': 8,
            }
        }

        with patch.object(InfoParser, 'data', new_callable=PropertyMock) as mock_data:
            mock_data.return_value = source
            result = parser.get_k_mesh()

        assert result == [
            {
                'grid': [4, 4, 4],
                'n_total_points': 64,
                'n_reduced_points': 8,
                'dimensionality': 3,
            }
        ]
        print(result)

    def test_maps_is_spin_polarization(self):
        parser = InfoParser()
        assert parser.get_is_spin_polarized('spin-unpolarised') is False
        assert parser.get_is_spin_polarized('spin-polarised') is True


@pytest.mark.unit
class TestBandstructureXMLParserMapping:
    def test_maps_k_path(self):
        parser = BandstructureXMLParser()
        source = {
            'bandstructure': {
                'band': [
                    {
                        'point': [
                            {'@coord': '0 0 0', '@eval': '-1'},
                            {'@coord': '0.5 0 0', '@eval': '1'},
                        ]
                    }
                ],
                'vertex': [{'@label': 'G', '@coord': '0 0 0'}],
            }
        }

        path = parser.get_k_path(source)

        assert path['n_line_points'] == 2
        np.testing.assert_allclose(path['points'], [[0, 0, 0], [0.5, 0, 0]])
        assert path['high_symmetry_path_names'] == ['G']


@pytest.mark.unit
class TestDosXMLParserMapping:
    def test_maps_energy_axis_and_values(self):
        parser = DosXMLParser()

        dos = parser.get_dos(
            {'point': [{'@e': '-1', '@dos': '2'}, {'@e': '0', '@dos': '4'}]}
        )

        np.testing.assert_allclose(dos['energy'], [-1, 0])
        np.testing.assert_allclose(dos['dos'], [2, 4])


@pytest.mark.unit
class TestEigvalParserMapping:
    @pytest.mark.parametrize(
        ('eigs_occs', 'expected'),
        [
            (
                [
                    {
                        'eigenvalues': [[-5.0, -3.0, 2.0]],
                        'occupancies': [[1.0, 1.0, 0.0]],
                    },
                    {
                        'eigenvalues': [[-4.5, -3.5, 1.0]],
                        'occupancies': [[1.0, 1.0, 0.0]],
                    },
                ],
                [dict(value=4.0)],
            ),
            (
                [
                    {
                        'eigenvalues': [[-5.0, -3.0, 2.0], [-4.0, -2.0, 3.0]],
                        'occupancies': [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
                    }
                ],
                [dict(value=5.0, spin_channel=0), dict(value=5.0, spin_channel=1)],
            ),
            ([{'eigenvalues': [[-5.0, -3.0]], 'occupancies': [[1.0, 1.0]]}], []),
            ([], []),
        ],
    )
    def test_maps_eigenvalues_to_band_gaps(self, eigs_occs, expected):
        parser = EigvalParser()

        gaps = parser.get_band_gaps({'eigenvalues_occupancies': eigs_occs})

        assert len(gaps) == len(expected)
        for gap, reference in zip(gaps, expected):
            assert gap['value'] == pytest.approx(reference['value'])
            assert gap.get('spin_channel') == reference.get('spin_channel')
