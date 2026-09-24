import numpy as np
import pytest

from nomad_simulation_parsers.parsers.vasp.chgcar_parser import CHGCARParser
from nomad_simulation_parsers.parsers.vasp.common import functional_key_from_params
from nomad_simulation_parsers.parsers.vasp.doscar_parser import DOSCARParser
from nomad_simulation_parsers.parsers.vasp.outcar_parser import OutcarParser
from nomad_simulation_parsers.parsers.vasp.xml_parser import VasprunParser
from tests.parsers.common import assert_approx


@pytest.mark.unit
class TestOutcarMapping:
    @pytest.mark.parametrize(
        ('parameters', 'expected'),
        [
            ({'GGA': 'PE'}, 'PBE'),
            ({'GGA': '--'}, 'PBE'),
            ({'METAGGA': 'SCAN'}, 'SCAN'),
            ({'LEXCH': 'CA'}, 'PZ81'),
            ({'LHFCALC': True, 'HFSCREEN': 0.2}, 'HSE06'),
            ({'LHFCALC': True, 'HFSCREEN': 0.3}, 'HSE03'),
            ({'LHFCALC': True, 'HFSCREEN': 0.2, 'GGA': 'PS'}, 'HSEsol'),
            ({'LHFCALC': True, 'HFSCREEN': 0.2, 'GGA': 'RP'}, None),
            ({'LHFCALC': True, 'AEXX': 0.25}, 'PBE0'),
            ({'LHFCALC': True, 'GGA': 'PBE'}, 'PBE0'),
            ({'LHFCALC': True, 'AEXX': 0.25, 'HFSCREEN': 0.5}, None),
            ({'LHFCALC': True, 'GGA': 'B5'}, 'B3LYP5'),
            ({'LHFCALC': True, 'GGA': 'B3'}, 'B3LYP'),
            (
                {'LHFCALC': True, 'AEXX': 1.0, 'ALDAC': 0.0, 'AGGAC': 0.0},
                None,
            ),
        ],
    )
    def test_maps_xc_tags_to_canonical_functionals(self, parameters, expected):
        assert functional_key_from_params(parameters) == expected

    def test_maps_scf_steps_with_energy_units(self):
        parser = OutcarParser()
        source = {
            'scf_iteration': [
                {'energy_total': -10.0, 'time': [1.0, 2.0]},
                {'energy_total': -10.5, 'time': [3.0, 4.0]},
            ]
        }

        result = parser.get_scf_steps(source)

        assert len(result['energies_total']) == 2
        assert result['energies_total'][-1].to('eV').magnitude == pytest.approx(-10.5)
        assert result['delta_energies_total'][-1].to('eV').magnitude == pytest.approx(
            0.5
        )
        assert_approx(result['durations'], [2.0, 4.0])

    def test_maps_common_source_helpers(self):
        parser = OutcarParser()
        forces = np.arange(6, dtype=float).reshape(2, 3)
        source = {'positions_forces': [np.zeros((2, 3)), forces]}

        assert (
            parser.get_version(
                {'version': '6.4.2', 'subversion': 'complex', 'platform': 'Linux'}
            )
            == '6.4.2 complex Linux'
        )
        assert parser.get_data({'value': 3.5}) == 3.5
        mapped_forces = parser.get_forces(source)
        assert_approx(mapped_forces['forces'], forces)
        assert mapped_forces['npoints'] == 2
        assert mapped_forces['rank'] == [3]
        assert parser.get_energy_contributions({'a': 1, 'b': 2}, exclude=['b']) == [
            {'name': 'a', 'value': 1}
        ]
        assert parser.get_functional_key({'GGA': 'PE'}) == 'PBE'
        assert parser.get_configurations({'step': 1})[0]['frame_index'] == 0
        assert parser.get_atoms([2, 1], [['PAW', 'Si'], ['PAW', 'O']]) == [
            {'label': 'Si'},
            {'label': 'Si'},
            {'label': 'O'},
        ]
        assert parser.get_atoms([1], [['PAW', 'Si']], frame_index=1) == []
        assert parser.get_periodic_boundary_conditions() == [True, True, True]

    def test_maps_eigenvalues_and_band_gaps(self):
        parser = OutcarParser()
        eigenvalues = np.array(
            [
                [0.0, -1.0, 1.0, 1.0, 0.5, 0.0],
                [0.0, -0.8, 1.0, 1.0, 0.7, 0.0],
            ]
        )

        mapped = parser.get_eigenvalues(eigenvalues, {'ISPIN': 2})
        gaps = parser.get_band_gaps(eigenvalues, {'ISPIN': 2})

        assert len(mapped) == 2
        assert [entry['spin_channel'] for entry in mapped] == [0, 1]
        assert len(gaps) == 2


@pytest.mark.unit
class TestDOSCARMapping:
    def test_maps_total_dos_to_spin_sections(self, tmp_path):
        doscar = tmp_path / 'DOSCAR'
        doscar.write_text(
            'header 1\n'
            'header 2\n'
            'header 3\n'
            'header 4\n'
            'header 5\n'
            '0.0 0.0 3 0.5 0.0\n'
            '-1.0 1.0 0.0\n'
            '0.0 2.0 0.0\n'
            '1.0 3.0 0.0\n'
        )
        outcar = tmp_path / 'OUTCAR'
        outcar.write_text('mock')
        parser = DOSCARParser()
        parser.filepath = str(outcar)

        dos = parser.get_dos(parser.data['total_dos'], parser.data['projected_dos'])

        assert len(dos) == 1
        assert dos[0]['spin'] is None
        assert_approx(dos[0]['dos'], [1.0, 2.0, 3.0])

    @pytest.mark.parametrize(
        ('outcar_name', 'doscar_name'),
        [
            ('OUTCAR', 'DOSCAR'),
            ('OUTCAR.relax', 'DOSCAR.relax'),
            ('OUTCAR.TR', 'DOSCAR.TR'),
        ],
    )
    def test_resolves_matching_doscar_suffix(self, tmp_path, outcar_name, doscar_name):
        doscar_lines = [
            'header0',
            'header1',
            'header2',
            'header3',
            'header4',
            '0.0 0.0 3 0.5 0.0',
            '-1.0 1.0 0.0',
            '0.0 2.0 0.0',
            '1.0 3.0 0.0',
        ]
        (tmp_path / doscar_name).write_text('\n'.join(doscar_lines))
        if doscar_name != 'DOSCAR':
            decoy_lines = doscar_lines[:6] + [
                '-1.0 9.0 0.0',
                '0.0 9.0 0.0',
                '1.0 9.0 0.0',
            ]
            (tmp_path / 'DOSCAR').write_text('\n'.join(decoy_lines))
        outcar = tmp_path / outcar_name
        outcar.write_text('dummy\n')

        parser = DOSCARParser()
        parser.filepath = str(outcar)
        dos = parser.get_dos(parser.data['total_dos'], parser.data['projected_dos'])

        assert len(dos) == 1
        assert parser.data['e_fermi'] == pytest.approx(0.5)
        assert list(dos[0]['dos']) == pytest.approx([1.0, 2.0, 3.0])


@pytest.mark.unit
class TestCHGCARMapping:
    def test_maps_auxiliary_charge_density(self, tmp_path):
        (tmp_path / 'OUTCAR').write_text('mock')
        (tmp_path / 'CHGCAR').write_text(
            'CHGCAR mock\n'
            '1.0\n'
            '1 0 0\n'
            '0 1 0\n'
            '0 0 1\n'
            '1\n'
            'Direct\n'
            '0 0 0\n\n'
            '2 1 1\n'
            '1.0 2.0\n'
        )
        parser = CHGCARParser()
        parser.filepath = str(tmp_path / 'OUTCAR')

        source = parser.to_dict()

        assert len(source['values']) == 1
        assert_approx(source['values'][0], [[[1.0]], [[2.0]]])


@pytest.mark.unit
class TestVasprunMapping:
    def test_maps_spin_sets_to_eigenvalues_and_occupations(self):
        parser = VasprunParser()
        source = {
            'set': [
                {'set': [{'r': [[-1.0, 1.0], [0.5, 0.0]]}]},
                {'set': [{'r': [[-0.8, 1.0], [0.7, 0.0]]}]},
            ]
        }

        result = parser.get_eigenvalues(source)

        assert len(result) == 2
        assert_approx(result[0]['value'], [[-1.0, 0.5]])
        assert_approx(result[1]['occupation'], [[1.0, 0.0]])
        assert [entry['spin_channel'] for entry in result] == [0, 1]

    def test_maps_dos_and_band_gaps(self):
        parser = VasprunParser()
        dos_source = {
            'i': {'@name': 'efermi', '__value': '1.5'},
            'total': {'array': {'set': [{'r': [[0.0, -2.0], [1.0, 3.0]]}]}},
        }
        eigenvalue_source = {
            'set': [{'r': [[-1.0, 1.0], [0.5, 0.0]]}],
        }

        dos = parser.get_total_dos(dos_source)
        gaps = parser.get_band_gaps(eigenvalue_source)

        assert_approx(dos[0]['energies'], [0.0, 1.0])
        assert_approx(dos[0]['value'], [2.0, 3.0])
        assert dos[0]['energy_fermi'] == 1.5
        assert gaps[0]['value'] == pytest.approx(1.5)

    def test_maps_common_source_helpers_and_parameters(self):
        parser = VasprunParser()
        parser._data = {
            'modeling': {
                'parameters': {
                    'separator': [
                        {
                            '@name': 'electronic',
                            'i': [
                                {'@name': 'GGA', '__value': 'PE'},
                                {'@name': 'EDIFF', '__value': '0.0001'},
                            ],
                        }
                    ]
                }
            }
        }

        assert parser.mix_alpha(0.5, True) == 0.5
        assert parser.mix_alpha(0.5, False) == 0
        assert parser.get_data({'__value': 2.0}) == 2.0
        assert parser.reshape_array(np.arange(6), (3,)).shape == (2, 3)
        assert parser._get_parameter('GGA', 'electronic') == 'PE'
        assert parser._find_parameters(('GGA', 'EDIFF')) == {
            'GGA': 'PE',
            'EDIFF': 0.0001,
        }
        assert parser.get_functional_key() == 'PBE'
        assert parser.get_configurations({'step': 1})[0]['frame_index'] == 0
        assert parser.get_atoms(
            [{'@name': 'atoms', 'set': {'rc': [{'c': ['Si']}, {'c': ['O']}]}}]
        ) == [{'label': 'Si'}, {'label': 'O'}]
        assert parser.get_atoms([], frame_index=1) == []
        assert parser.get_positions([[0.5, 0.0, 0.0]], np.eye(3)).tolist() == [
            [0.5, 0.0, 0.0]
        ]
        assert parser.get_positions(None) is None
        assert parser.get_periodic_boundary_conditions() == [True, True, True]

    def test_maps_scf_steps_and_energy_contributions(self):
        parser = VasprunParser()
        source = {
            'scstep': [
                {
                    'energy': {
                        'i': [
                            {'@name': 'e_fr_energy', '__value': '-10.0'},
                        ]
                    },
                    'time': [{'@name': 'total', '__value': [1.0, 2.0]}],
                },
                {
                    'energy': {
                        'i': [
                            {'@name': 'e_fr_energy', '__value': '-10.5'},
                        ]
                    },
                    'time': [{'@name': 'total', '__value': [3.0, 4.0]}],
                },
            ]
        }

        result = parser.get_scf_steps(source)

        assert len(result['energies_total']) == 2
        assert result['energies_total'][-1].to('eV').magnitude == pytest.approx(-10.5)
        assert result['delta_energies_total'][-1].to('eV').magnitude == pytest.approx(
            0.5
        )
        assert result['durations'] == [2.0, 4.0]
        assert parser.get_energy_contributions(
            [{'@name': 'a'}, {'@name': 'b'}], exclude=['b']
        ) == [{'@name': 'a'}]
