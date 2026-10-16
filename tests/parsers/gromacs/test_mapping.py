from pathlib import Path

import numpy as np
import pytest

from nomad_simulation_parsers.parsers.gromacs.mdp_parser import (
    GromacsMdpParser as GromacsMdpFileParser,
)
from nomad_simulation_parsers.parsers.gromacs.parser import (
    GromacsEDRParser,
    GromacsLogParser,
    GromacsMDAnalysisParser,
    GromacsMDPParser,
    GromacsMetainfoParser,
    GromacsXVGParser,
    Program,
    Simulation,
)
from nomad_simulation_parsers.parsers.gromacs.xvg_parser import (
    GromacsXvgParser as GromacsXvgFileParser,
)
from tests.parsers.common import assert_approx


class StubTrajectory:
    def __init__(self):
        self.positions = [np.zeros((2, 3)), np.ones((2, 3))]

    def get_positions(self, index):
        return self.positions[index]

    def get_velocities(self, index):
        return np.full((2, 3), index, dtype=float)

    def get_lattice_vectors(self, index):
        return np.eye(3) * (index + 1)

    def get_n_atoms(self, index):
        return self.positions[index].shape[0]

    def get_atom_labels(self, index):
        return ['O', 'H']

    def get_frame_data(self, index):
        return {
            'positions': self.get_positions(index),
            'velocities': self.get_velocities(index),
            'lattice_vectors': self.get_lattice_vectors(index),
        }

    def get_interactions(self):
        return [{'type': 'O-H', 'atom_indices': [0, 1], 'atom_labels': ['O', 'H']}]

    def get(self, key, default=None):
        return default


class StubInteractionsDataObject:
    def __init__(self, interactions):
        self._interactions = interactions

    def get_interactions(self):
        return self._interactions

    def get(self, key, default=None):
        return default


class StubMDAnalysisDataObject:
    def __init__(self, positions, velocities, lattices):
        self._positions = positions
        self._velocities = velocities
        self._lattices = lattices

    def get_positions(self, index):
        return np.asarray(self._positions[index])

    def get_velocities(self, index):
        return None if self._velocities is None else np.asarray(self._velocities[index])

    def get_lattice_vectors(self, index):
        return None if self._lattices is None else np.asarray(self._lattices[index])

    def get_n_atoms(self, index):
        return self._positions[index].shape[0]

    def get_atom_labels(self, index):
        return ['H'] * self.get_n_atoms(index)

    def get_frame_data(self, index):
        return {
            'positions': self.get_positions(index),
            'velocities': self.get_velocities(index),
            'lattice_vectors': self.get_lattice_vectors(index),
        }

    def get_interactions(self):
        return []

    def get(self, key, default=None):
        return default


@pytest.mark.unit
class TestGromacsMDAnalysisMapping:
    @pytest.fixture
    def parser(self):
        return GromacsMDAnalysisParser()

    def test_maps_trajectory_frames_and_topology_once(self, parser):
        parser.data_object = StubTrajectory()
        parser._trajectory_steps_sampled = [0, 1]

        configurations = parser.get_configurations()

        assert len(configurations) == 2
        assert [label['label'] for label in configurations[0]['labels']] == ['O', 'H']
        assert 'labels' not in configurations[1]
        assert configurations[0]['positions'].shape == (2, 3)
        assert_approx(configurations[0]['bond_list'], [[0, 1]])

    def test_maps_force_field_contributions_by_interaction_type(self, parser):
        parser.data_object = StubTrajectory()

        contributions = parser.get_force_field_contributions()

        assert contributions == [
            {
                'functional_form': 'O-H',
                'particle_indices': [[0, 1]],
                'particle_labels': [['O', 'H']],
            }
        ]

    def test_maps_configurations_from_stub_data(self, parser):
        parser._trajectory_steps_sampled = [0, 1]
        parser.data_object = StubMDAnalysisDataObject(
            [np.zeros((3, 3)), np.ones((3, 3))],
            [np.zeros((3, 3)), np.ones((3, 3)) * 2],
            [np.eye(3), np.eye(3) * 2],
        )

        configurations = parser.get_configurations()

        assert len(configurations) == 2
        assert_approx(configurations[0]['positions'], np.zeros((3, 3)))
        assert_approx(configurations[1]['positions'], np.ones((3, 3)))
        assert_approx(configurations[0]['lattice_vectors'], np.eye(3))
        assert 'labels' in configurations[0]

    def test_converts_stub_configurations_to_simulation(self, parser):
        parser._trajectory_steps_sampled = [0, 1]
        parser.data_object = StubMDAnalysisDataObject(
            [np.zeros((2, 3)), np.ones((2, 3))],
            [None, None],
            [np.eye(3), np.eye(3)],
        )
        simulation = Simulation(program=Program(name='GROMACS'))
        metainfo_parser = GromacsMetainfoParser()
        metainfo_parser.data_object = simulation
        metainfo_parser.annotation_key = 'TPR'
        parser.convert(metainfo_parser)

        if simulation.model_system:
            assert len(simulation.model_system) == 2
            assert simulation.model_system[0].positions is not None
        else:
            assert len(parser.get_configurations()) == 2

    def test_groups_force_field_interactions(self, parser):
        parser.data_object = StubInteractionsDataObject(
            [
                {
                    'type': 'bond_harmonic',
                    'atom_indices': [0, 1],
                    'atom_labels': ['O', 'H'],
                },
                {
                    'type': 'bond_harmonic',
                    'atom_indices': [2, 3],
                    'atom_labels': ['O', 'H'],
                },
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
        )

        contributions = parser.get_force_field_contributions()

        assert {item['functional_form'] for item in contributions} == {
            'bond_harmonic',
            'angle_harmonic',
        }
        assert all(item['particle_indices'] for item in contributions)

    @pytest.mark.parametrize(
        'interactions,expected',
        [
            ([{'atom_indices': [0, 1]}], (1, 2)),
            ([{'atom_indices': [0, 1, 2]}], None),
            ([], None),
            ([{'atom_indices': [0, 1]}, {'atom_indices': None}], (1, 2)),
        ],
    )
    def test_extracts_valid_bonds_only(self, parser, interactions, expected):
        parser.data_object = StubInteractionsDataObject(interactions)
        bonds = parser.get_bond_list()

        if expected is None:
            assert bonds is None
        else:
            assert bonds.shape == expected


@pytest.mark.unit
class TestGromacsLogMapping:
    @pytest.fixture
    def parser(self):
        return GromacsLogParser()

    @pytest.mark.parametrize(
        ('integrator', 'expected'),
        [
            ('md', 'leap_frog'),
            ('md-vv', 'velocity_verlet'),
            ('steep', 'steepest_descent'),
        ],
    )
    def test_maps_gromacs_integrator(self, parser, integrator, expected):
        assert parser.get_integrator_type(integrator) == expected

    def test_maps_thermodynamic_outputs(self, parser):
        parser._data = {
            'Time': [0.0],
            'Potential': [-10.0],
            'Total Energy': [-8.0],
            'Temperature': [300.0],
        }
        parser._thermodynamic_steps = [0]
        parser._trajectory_steps = [0]

        output = parser.get_outputs()[0]

        assert output['step'] == 0
        assert output['system_ref'] == '/data/model_system/0'
        assert output['temperatures'][0]['name'] == 'Temperature'
        assert_approx(output['energy']['value'].magnitude, -8.0)

    def test_maps_periodic_boundary_configurations(self, parser):
        parser._data = {'input_parameters': {'pbc': 'xy'}}
        parser._trajectory_steps_sampled = [0, 1, 2]
        configurations = parser.get_configurations()
        assert len(configurations) == 3
        assert all(config['pbc'] == [True, True, False] for config in configurations)

    def test_maps_coulomb_types(self, parser):
        values = {
            'cut-off': 'cutoff',
            'cutoff': 'cutoff',
            'Ewald': 'ewald',
            'PME': 'particle_mesh_ewald',
            'P3M-AD': 'particle_particle_particle_mesh',
            'Reaction-Field': 'reaction_field',
            'Reaction-Field-zero': 'reaction_field',
            'unknown': None,
        }
        for source, expected in values.items():
            assert parser.get_coulomb_type(source) == expected

    def test_prefers_compressed_coordinate_frequency(self, parser):
        assert (
            parser.get_coordinate_save_frequency(
                {'nstxout-compressed': 100, 'nstxout': 50}
            )
            == 100
        )
        assert parser.get_coordinate_save_frequency({'nstxout': 50}) == 50
        assert (
            parser.get_coordinate_save_frequency(
                {'nstxout-compressed': 0, 'nstxout': 50}
            )
            == 50
        )
        assert parser.get_coordinate_save_frequency({}) is None
        assert parser.get_coordinate_save_frequency(None) is None

    @pytest.mark.parametrize(
        'params,expected',
        [
            ({'tcoupl': 'v-rescale', 'pcoupl': 'parrinello-rahman'}, 'NPT'),
            ({'tcoupl': 'nose-hoover', 'pcoupl': 'no'}, 'NVT'),
            ({'tcoupl': 'no', 'pcoupl': 'berendsen'}, 'NPH'),
            ({'tcoupl': 'no', 'pcoupl': 'no'}, 'NVE'),
            ({}, 'NVE'),
            (None, None),
        ],
    )
    def test_maps_thermodynamic_ensemble(self, parser, params, expected):
        assert parser.get_thermodynamic_ensemble(params) == expected

    @pytest.mark.parametrize(
        'source,expected',
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
    def test_maps_thermostat_type(self, parser, source, expected):
        assert parser.get_thermostat_type(source) == expected

    def test_reads_thermostat_parameters(self, parser):
        assert parser.get_reference_temperature({'grpopts': {'ref-t': 300.0}}) == 300.0
        assert (
            parser.get_reference_temperature({'grpopts': {'ref-t': [300.0, 310.0]}})
            == 300.0
        )
        assert parser.get_reference_temperature({'grpopts': {'ref-t': []}}) is None
        assert (
            parser.get_thermostat_coupling_constant({'grpopts': {'tau-t': [0.1, 0.2]}})
            == 0.1
        )
        assert parser.get_thermostat_coupling_constant({}) is None

    @pytest.mark.parametrize(
        'source,expected',
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
    def test_maps_barostat_types(self, parser, source, expected):
        assert parser.get_barostat_type(source) == expected

    @pytest.mark.parametrize(
        'source,expected',
        [
            ('isotropic', 'isotropic'),
            ('semiisotropic', 'semi_isotropic'),
            ('anisotropic', 'anisotropic'),
            ('surface-tension', 'surface_tension'),
            ('unknown_type', None),
            (None, None),
        ],
    )
    def test_maps_barostat_coupling_types(self, parser, source, expected):
        assert parser.get_barostat_coupling_type(source) == expected

    def test_reads_barostat_parameters_and_matrix(self, parser):
        assert parser.get_barostat_coupling_constant({'tau-p': 2.0}) == 2.0
        assert parser.get_barostat_coupling_constant({}) is None
        matrix = np.eye(3)
        assert_approx(parser.get_matrix_parameter({'ref-p': matrix}, 'ref-p'), matrix)
        assert parser.get_matrix_parameter({'ref-p': np.eye(2)}, 'ref-p') is None
        assert parser.get_matrix_parameter(None, 'ref-p') is None

    def test_maps_free_energy_parameters(self, parser):
        assert parser.get_free_energy_calc_type({'free-energy': 'yes'}) == 'alchemical'
        assert (
            parser.get_free_energy_calc_type({'free-energy': 'umbrella'})
            == 'umbrella_sampling'
        )
        assert parser.get_free_energy_calc_type({'free-energy': 'no'}) is None
        source = {'input_parameters': {'free-energy': 'yes'}}
        assert parser.get_fep_params_if_active(source) is source
        assert parser.get_fep_params_if_active({'input_parameters': {}}) is None
        assert parser.get_lambda_state_index({'init-lambda-state': '3'}) == 3
        assert parser.get_lambda_state_index({'init-lambda-state': -1}) is None

    def test_maps_lambda_schedules(self, parser):
        params = {
            'all-lambdas': {
                'vdw-lambdas': '0.0 0.25 0.5 0.75 1.0',
                'coul-lambdas': '0.0 0.0 0.0 0.0 0.0',
            },
            'sc-alpha': '0.5',
            'sc-power': '1',
            'sc-sigma': '0.3',
        }
        result = parser.get_lambdas_schedule(params)
        assert len(result) == 1
        assert result[0]['interaction_type'] == 'vdw'
        assert_approx(result[0]['lambda_values'], [0.0, 0.25, 0.5, 0.75, 1.0])
        assert result[0]['softcore_enabled'] is True
        assert_approx(
            parser.get_current_lambdas(
                {
                    'init-lambda-state': '2',
                    'all-lambdas': {'vdw-lambdas': '0.0 0.5 1.0'},
                    'sc-alpha': '0.0',
                }
            ),
            [1.0],
        )


@pytest.mark.unit
class TestGromacsMDPMapping:
    @pytest.fixture
    def parser(self):
        return GromacsMDPParser(text_parser=GromacsMdpFileParser())

    def test_adds_unique_mdp_parameters(self, tmp_path, parser):
        mdp = tmp_path / 'md.mdp'
        mdp.write_text('integrator = md\ncustom-setting = BAR\n')
        parser.filepath = str(mdp)
        source = parser.load_file()
        source.parse()

        parameters = source.results['input_parameters']
        assert parameters['custom-setting'] == 'BAR'
        assert 'mdp_unique_params' in parameters


@pytest.mark.unit
class TestGromacsEDRMapping:
    @pytest.fixture
    def parser(self):
        return GromacsEDRParser()

    def test_maps_edr_energies(self, parser):
        parser._data = {
            'Time': [0.0, 1.0],
            'Potential': [10.0, 20.0],
            'Kinetic En.': [5.0, 6.0],
            'Total Energy': [15.0, 26.0],
        }
        parser._thermodynamic_steps = [0, 1]
        energies = parser.get_energies()
        assert len(energies) == 2
        assert all(hasattr(energy, 'units') for energy in energies)

    def test_maps_energies_to_outputs(self, parser):
        parser._data = {
            'Time': [0.0, 1.0],
            'Potential': [10.0, 20.0],
            'Kinetic En.': [5.0, 6.0],
            'Total Energy': [15.0, 26.0],
        }
        parser._thermodynamic_steps = [0, 1]

        energies = parser.get_energies()
        outputs = parser.get_outputs()

        assert len(energies) == 2
        assert len(outputs) == 2
        assert all(hasattr(value, 'units') for value in energies)
        assert_approx(energies[0].magnitude, 15.0)


@pytest.mark.unit
class TestGromacsXVGMapping:
    @pytest.fixture
    def parser(self):
        return GromacsXVGParser(text_parser=GromacsXvgFileParser())

    def test_maps_free_energy_columns(self, parser):
        parser._data = {
            'title': r'dH/d\xl\f{} and \xD\f{}H',
            'xaxis': 'Time (ps)',
            'column_vals': np.array(
                [
                    [0.0, 10.0, 1.0, 2.0, 3.0, 4.0],
                    [1.0, 11.0, 1.5, 2.5, 3.5, 4.5],
                ]
            ),
        }

        results = parser.get_results()['free_energy_calculations']

        assert results['n_frames'] == 2
        assert results['n_states'] == 2
        assert_approx(results['value_total_energy'].magnitude, [10.0, 11.0])
        assert_approx(parser.get_fep_xvg_data('n_frames'), 2)

    def test_xvg_parser_get_results_valid(self, parser):
        path = (
            Path(__file__).resolve().parents[2]
            / 'data'
            / 'gromacs'
            / 'free_energy_calculations'
            / 'alchemical_transformation_single_run'
            / 'fep_run-7.xvg'
        )
        parser.text_parser.mainfile = str(path)
        parser.text_parser.parse()
        parser.filepath = str(path)
        free_energy = parser.get_results()['free_energy_calculations']
        assert free_energy['n_frames'] == 5001
        assert free_energy['n_states'] == 11
        assert free_energy['times'] is not None
        assert free_energy['value_total_energy_derivative'].shape == (5001,)
        assert free_energy['value_total_energy_differences'].shape == (5001, 11)
        assert free_energy['value_PV_energy'].shape == (5001,)
