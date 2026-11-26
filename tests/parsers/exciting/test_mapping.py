from unittest.mock import PropertyMock, patch

import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.exciting.parser import (
    BandstructureXMLParser,
    DosXMLParser,
    EigvalParser,
    InfoParser,
    SpeciesXMLParser,
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


@pytest.mark.unit
class TestSpeciesXMLParserMapping:
    @pytest.mark.parametrize(
        ('basis', 'lmaxapw', 'expected'),
        [
            (None, 2, []),
            ({}, 2, []),
            ({'default': {'@type': 'apw'}}, None, []),
            (
                {'default': {'@type': 'lapw'}},
                1,
                [
                    {
                        'name': 0,
                        'orbitals': [{'@type': 'lapw', '@l': '0'}],
                    },
                    {
                        'name': 1,
                        'orbitals': [{'@type': 'lapw', '@l': '1'}],
                    },
                ],
            ),
            (
                {
                    'default': {'@type': 'lapw'},
                    'custom': [{'@l': '1', '@type': 'apw'}],
                    'lo': {
                        '@l': '1',
                        'wf': {'@matchingOrder': '0'},
                    },
                },
                1,
                [
                    {
                        'name': 0,
                        'orbitals': [{'@type': 'lapw', '@l': '0'}],
                    },
                    {
                        'name': 1,
                        'orbitals': [
                            {'@l': '1', '@type': 'apw'},
                            {
                                '@l': '1',
                                'wf': {'@matchingOrder': '0'},
                            },
                        ],
                    },
                ],
            ),
            (
                {'default': {'@type': 'apw+lo'}},
                0,
                [
                    {
                        'name': 0,
                        'orbitals': [
                            {'@type': 'apw', '@l': '0'},
                            {
                                '@l': '0',
                                'wf': [
                                    {'@matchingOrder': '0'},
                                    {'@matchingOrder': '1'},
                                ],
                            },
                        ],
                    }
                ],
            ),
        ],
        ids=[
            'missing-basis',
            'empty-basis',
            'missing-lmaxapw',
            'default-basis',
            'custom-basis-and-local-orbital',
            'apw-plus-lo',
        ],
    )
    def test_builds_l_channels(self, basis, lmaxapw, expected):
        parser = SpeciesXMLParser()
        parser.lmaxapw = lmaxapw

        assert parser.build_l_channels(basis) == expected

    @pytest.mark.parametrize(
        ('apw', 'l_value', 'expected'),
        [
            (
                {
                    '@type': 'apw+lo',
                    '@trialEnergy': '0.15',
                    '@searchE': 'true',
                },
                2,
                [
                    {
                        '@type': 'apw',
                        '@l': '2',
                        '@trialEnergy': '0.15',
                        '@searchE': 'true',
                    },
                    {
                        '@l': '2',
                        'wf': [
                            {
                                '@matchingOrder': '0',
                                '@trialEnergy': '0.15',
                                '@searchE': 'true',
                            },
                            {
                                '@matchingOrder': '1',
                                '@trialEnergy': '0.15',
                                '@searchE': 'true',
                            },
                        ],
                    },
                ],
            ),
            (
                {'@type': 'apw+lo'},
                0,
                [
                    {'@type': 'apw', '@l': '0'},
                    {
                        '@l': '0',
                        'wf': [
                            {'@matchingOrder': '0'},
                            {'@matchingOrder': '1'},
                        ],
                    },
                ],
            ),
        ],
        ids=['with-energy-parameters', 'without-energy-parameters'],
    )
    def test_expands_apw_plus_lo(self, apw, l_value, expected):
        parser = SpeciesXMLParser()
        assert parser._expand_apw_plus_lo(apw, l_value) == expected

    @pytest.mark.parametrize(
        ('orbital', 'expected'),
        [
            (
                {'@type': 'apw', '@trialEnergy': '0.15'},
                [0.15],
            ),
            (
                {'@type': 'lapw', '@trialEnergy': '0.15'},
                [0.15, 0.15],
            ),
            (
                {'wf': {'@trialEnergy': '-0.25'}},
                [-0.25],
            ),
            (
                {
                    'wf': [
                        {'@trialEnergy': '-0.25'},
                        {'@trialEnergy': '0.40'},
                    ]
                },
                [-0.25, 0.40],
            ),
            (
                {
                    'wf': [
                        {'@trialEnergy': '0.15'},
                        {'@matchingOrder': '1'},
                    ]
                },
                None,
            ),
            ({'@type': 'apw'}, None),
            (None, None),
        ],
        ids=[
            'apw',
            'lapw',
            'single-wavefunction',
            'multiple-wavefunctions',
            'incomplete-wavefunctions',
            'missing-trial-energy',
            'absent',
        ],
    )
    def test_maps_orbital_energy_parameter(self, orbital, expected):
        parser = SpeciesXMLParser()

        result = parser.get_orbital_energy_parameter(orbital)

        if expected is None:
            assert result is None
        else:
            np.testing.assert_allclose(
                result.to('hartree').magnitude,
                expected,
            )

    @pytest.mark.parametrize(
        ('orbital', 'expected'),
        [
            (
                {'@type': 'apw', '@searchE': 'false'},
                ['fixed'],
            ),
            (
                {'@type': 'lapw', '@searchE': 'true'},
                ['initial', 'initial'],
            ),
            (
                {'wf': {'@searchE': 'false'}},
                ['fixed'],
            ),
            (
                {
                    'wf': [
                        {'@searchE': 'false'},
                        {'@searchE': 'true'},
                    ]
                },
                ['fixed', 'initial'],
            ),
            (
                {'@type': 'apw', '@searchE': 'unknown'},
                None,
            ),
            ({'@type': 'lapw'}, None),
            (None, None),
        ],
        ids=[
            'apw-fixed',
            'lapw-initial',
            'single-wavefunction',
            'mixed-wavefunctions',
            'unknown-status',
            'missing-status',
            'absent',
        ],
    )
    def test_maps_energy_status(self, orbital, expected):
        parser = SpeciesXMLParser()
        assert parser.get_energy_status(orbital) == expected

    @pytest.mark.parametrize(
        ('orbital', 'expected'),
        [
            (
                {'@type': 'apw'},
                [0],
            ),
            (
                {'@type': 'lapw'},
                [0, 1],
            ),
            (
                {'wf': {'@matchingOrder': '2'}},
                [2],
            ),
            (
                {
                    'wf': [
                        {'@matchingOrder': '0'},
                        {'@matchingOrder': '1'},
                    ]
                },
                [0, 1],
            ),
            (
                {
                    'wf': [
                        {'@matchingOrder': '0'},
                        {},
                    ]
                },
                None,
            ),
            ({'@type': 'unknown'}, None),
            (None, None),
        ],
        ids=[
            'apw',
            'lapw',
            'single-wavefunction',
            'multiple-wavefunctions',
            'incomplete-wavefunctions',
            'unknown-orbital',
            'absent',
        ],
    )
    def test_maps_orbital_differential_order(self, orbital, expected):
        parser = SpeciesXMLParser()

        result = parser.get_orbital_differential_order(orbital)

        if expected is None:
            assert result is None
        else:
            np.testing.assert_array_equal(result, expected)
            assert result.dtype == np.int32

    @pytest.mark.parametrize(
        ('species_file', 'atoms', 'expected'),
        [
            (
                'C.xml',
                [
                    {'file': 'C.xml'},
                    {'file': 'O.xml'},
                    {'file': 'C.xml'},
                ],
                [
                    '#/data/model_system/0/particle_states/0',
                    '#/data/model_system/0/particle_states/2',
                ],
            ),
            (
                'I.xml',
                [
                    {'file': 'Pb.xml'},
                    {'file': 'I.xml'},
                    {'file': 'I.xml'},
                ],
                [
                    '#/data/model_system/0/particle_states/1',
                    '#/data/model_system/0/particle_states/2',
                ],
            ),
            (
                'C.xml',
                [{'file': 'O.xml'}],
                [],
            ),
            (
                'C.xml',
                [],
                [],
            ),
            (
                None,
                [{'file': 'C.xml'}],
                [],
            ),
            (
                'C.xml',
                None,
                [],
            ),
        ],
        ids=[
            'non-consecutive-atoms',
            'multiple-atoms-same-species',
            'no-matching-species',
            'empty-atoms',
            'absent-species-file',
            'absent-species-list',
        ],
    )
    def test_maps_muffin_tin_species_scope(self, species_file, atoms, expected):
        parser = InfoParser()

        species_list = (
            [{'file': atom['file']} for atom in atoms] if atoms is not None else None
        )

        with patch.object(
            parser, '_get_species_atoms', return_value=atoms
        ) as mock_atoms:
            result = parser.get_muffin_tin_species_scope(species_file, species_list)

            assert result == expected

            if species_file is None or species_list is None:
                mock_atoms.assert_not_called()
            else:
                mock_atoms.assert_called_once_with(species_list)

    @pytest.mark.parametrize(
        ('species_list', 'expected'),
        [
            (
                [
                    {
                        'name': 'C',
                        'file': 'C.xml',
                        'positions': [[0.0, 0.0, 0.0]],
                    },
                ],
                [
                    {'name': 'C', 'file': 'C.xml'},
                ],
            ),
            (
                [
                    {
                        'name': 'C',
                        'file': 'C.xml',
                        'positions': [
                            [0.0, 0.0, 0.0],
                            [0.5, 0.5, 0.5],
                        ],
                    },
                ],
                [
                    {'name': 'C', 'file': 'C.xml'},
                    {'name': 'C', 'file': 'C.xml'},
                ],
            ),
            (
                [
                    {
                        'name': 'Pb',
                        'file': 'Pb.xml',
                        'positions': [[0.0, 0.0, 0.0]],
                    },
                    {
                        'name': 'I',
                        'file': 'I.xml',
                        'positions': [
                            [0.5, 0.0, 0.0],
                            [0.0, 0.5, 0.0],
                        ],
                    },
                ],
                [
                    {'name': 'Pb', 'file': 'Pb.xml'},
                    {'name': 'I', 'file': 'I.xml'},
                    {'name': 'I', 'file': 'I.xml'},
                ],
            ),
            (
                [
                    {
                        'name': 'C',
                        'file': 'C_1.xml',
                        'positions': [[0.0, 0.0, 0.0]],
                    },
                    {
                        'name': 'C',
                        'file': 'C_2.xml',
                        'positions': [[0.5, 0.5, 0.5]],
                    },
                ],
                [
                    {'name': 'C', 'file': 'C_1.xml'},
                    {'name': 'C', 'file': 'C_2.xml'},
                ],
            ),
            (
                [
                    {
                        'name': 'C',
                        'file': 'C.xml',
                        'positions': [[0.0, 0.0, 0.0]],
                        'positions_format': 'lattice',
                        'radial_points': 300,
                        'muffin_tin_radius': 1.5,
                    },
                ],
                [
                    {
                        'name': 'C',
                        'file': 'C.xml',
                        'muffin_tin_radius': 1.5,
                    },
                ],
            ),
            (
                [
                    {'name': 'C', 'file': 'C.xml'},
                ],
                [],
            ),
            (
                [
                    {
                        'name': 'C',
                        'file': 'C.xml',
                        'positions': None,
                    },
                ],
                [],
            ),
            (
                [
                    {
                        'name': 'C',
                        'file': 'C.xml',
                        'positions': [],
                    },
                ],
                [],
            ),
            (
                [],
                [],
            ),
        ],
        ids=[
            'single-atom',
            'multiple-atoms',
            'multiple-species',
            'same-element-different-species-files',
            'exclude-species-fields',
            'missing-positions',
            'none-positions',
            'empty-positions',
            'empty-species-list',
        ],
    )
    def test_get_species_atoms(self, species_list, expected):
        assert InfoParser._get_species_atoms(species_list) == expected
