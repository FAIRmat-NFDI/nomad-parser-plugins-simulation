import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.exciting.parser import (
    BandstructureXMLParser,
    DosXMLParser,
    EigvalParser,
    InfoParser,
)


@pytest.mark.unit
class TestInfoParserMapping:
    @pytest.mark.parametrize(
        ('xc_type', 'expected_key', 'expected_components'),
        [
            (20, 'PBE', [{'libxc': 'GGA_C_PBE'}, {'libxc': 'GGA_X_PBE'}]),
            (
                22,
                'PBESOL',
                [{'libxc': 'GGA_C_PBE_SOL'}, {'libxc': 'GGA_X_PBE_SOL'}],
            ),
            (408, 'HSE06', [{'libxc': 'HYB_GGA_XC_HSE06'}]),
            (100, None, []),
            (0, None, []),
            (None, None, []),
        ],
        ids=['pbe', 'pbesol', 'hse06', 'libxc-from-input-xml', 'zero', 'absent'],
    )
    def test_maps_xc_type_to_functional(
        self, xc_type, expected_key, expected_components
    ):
        parser = InfoParser()
        assert parser.get_xc_functional_key(xc_type) == expected_key
        assert parser.get_xc_functional_components(xc_type) == expected_components

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
        ('smearing_kind', 'expected', 'warns'),
        [
            ('Gaussian', 'Gaussian', False),
            ('Fermi-Dirac', 'Fermi-Dirac', False),
            (
                ['Methfessel-Paxton', 'order', '1,', 'Phys.', 'Rev.', 'B', '40,']
                + [3616, '(1989)'],
                'Methfessel-Paxton 1',
                False,
            ),
            (
                ['Methfessel-Paxton', 'order', '2,', 'Phys.', 'Rev.', 'B', '40,']
                + [3616, '(1989)'],
                'Methfessel-Paxton 2',
                False,
            ),
            (['Extended', 'linear', 'tetrahedron', 'method'], 'Tetrahedra', False),
            (['Square-wave', 'impulse'], 'Square-wave', False),
            ('Unknown-scheme', None, True),
            ('', None, False),
            ([], None, False),
            (None, None, False),
        ],
        ids=[
            'gaussian',
            'fermi-dirac',
            'methfessel-paxton-order-1',
            'methfessel-paxton-order-2',
            'tetrahedra',
            'square-wave',
            'unknown-warns',
            'empty-string',
            'empty-tokens',
            'absent',
        ],
    )
    def test_maps_smearing_kind(self, caplog, smearing_kind, expected, warns):
        parser = InfoParser()
        assert parser.get_smearing_kind(smearing_kind) == expected
        assert [record.levelname for record in caplog.records] == (
            ['WARNING'] if warns else []
        )

    @pytest.mark.parametrize(
        ('spin_treatment', 'expected', 'warns'),
        [
            ('spin-polarised', True, False),
            ('spin-unpolarised', False, False),
            ('unknown', None, True),
            (None, None, False),
        ],
        ids=['polarised', 'unpolarised', 'unknown-warns', 'absent'],
    )
    def test_maps_spin_treatment(self, caplog, spin_treatment, expected, warns):
        parser = InfoParser()
        assert parser.get_is_spin_polarized(spin_treatment) is expected
        assert [record.levelname for record in caplog.records] == (
            ['WARNING'] if warns else []
        )

    @pytest.mark.parametrize(
        ('initialization', 'expected'),
        [
            (
                {'kpoint_grid': np.array([4, 4, 4]), 'x_exciting_number_kpoints': 8},
                [
                    {
                        'grid': [4, 4, 4],
                        'n_total_points': 64,
                        'n_reduced_points': 8,
                        'dimensionality': 3,
                    }
                ],
            ),
            (
                {'kpoint_grid': np.array([6, 6, 4])},
                [
                    {
                        'grid': [6, 6, 4],
                        'n_total_points': 144,
                        'n_reduced_points': None,
                        'dimensionality': 3,
                    }
                ],
            ),
            (
                {'kpoint_grid': np.array([1, 1, 1])},
                [
                    {
                        'grid': [1, 1, 1],
                        'n_total_points': 1,
                        'n_reduced_points': None,
                        'dimensionality': 3,
                    }
                ],
            ),
            ({'x_exciting_number_kpoints': 8}, []),
            ({}, []),
        ],
        ids=['cubic-grid', 'anisotropic-grid', 'gamma-only', 'no-grid', 'empty'],
    )
    def test_maps_k_mesh(self, initialization, expected):
        k_mesh = InfoParser().get_k_mesh(initialization)
        assert [{**mesh, 'grid': list(mesh['grid'])} for mesh in k_mesh] == expected


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
