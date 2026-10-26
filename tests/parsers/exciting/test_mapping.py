import numpy as np
import pytest
from unittest.mock import PropertyMock, patch
from nomad.units import ureg
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


    def test_maps_smearing_kind(self):
        parser = InfoParser()
        assert parser.get_smearing_kind('Extended linear tetrahedron method'
            ) == 'Tetrahedra-extended'

    def test_maps_k_mesh(self):
        parser = InfoParser()
        source = {
            'initialization': {
                'kpoint_grid': [4, 4, 4],
                'x_exciting_number_kpoints': 8,
            }
        }

        with patch.object(
            InfoParser, 'data', new_callable=PropertyMock
        ) as mock_data:
            mock_data.return_value = source
            result = parser.get_k_mesh()

        assert result  == [
            {'grid': [4, 4, 4], 
            'n_total_points': 64, 
            'n_reduced_points': 8, 
            'dimensionality': 3}
            ]
        print(result)

    def test_maps_is_spin_polarization(self):
        parser = InfoParser()
        assert not parser.get_is_spin_polarized('spin-unpolarised') 
        assert parser.get_is_spin_polarized('spin-polarised') 

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
