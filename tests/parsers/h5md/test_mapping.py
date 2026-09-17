import numpy as np
import pytest

from nomad_simulation_parsers.parsers.h5md.parser import H5MDH5Parser
from tests.parsers.common import approx, assert_approx


@pytest.mark.unit
class TestH5MDMapping:
    def test_maps_unit_factor(self):
        parser = H5MDH5Parser()
        source = {
            'value': {
                '__value': np.array([1.0, 2.0]),
                '@unit': 'angstrom',
                '@unit_factor': 2.0,
            }
        }

        value = parser.get_value('value', source)

        assert_approx(value.to('angstrom').magnitude, [2.0, 4.0])

    def test_gets_nested_source_data(self):
        parser = H5MDH5Parser(data={'particles': {'all': {'label': 'all'}}})

        assert parser.get_source(parser.data, 'particles.all.label') == 'all'
        assert parser.get_source(parser.data, None) == {}

    def test_maps_case_insensitive_values(self):
        parser = H5MDH5Parser()
        source = {'value': {'__value': 'npt'}}

        assert parser.map_value(source, 'value', 'upper') == 'NPT'
        assert parser.map_value(source, 'value', 'lower') == 'npt'

    def test_maps_subsystems_only_for_topology_frame(self):
        parser = H5MDH5Parser(
            data={
                'connectivity': {
                    'particles_group': {
                        'water': {'label': 'water'},
                        'ion': {'label': 'ion'},
                    }
                }
            }
        )

        assert parser.get_sub_systems({'step': 0}, path='connectivity') == [
            {'label': 'water'},
            {'label': 'ion'},
        ]
        assert parser.get_sub_systems({'step': 1}, path='connectivity') == []

    def test_maps_step_data_with_units(self):
        parser = H5MDH5Parser()
        data = {
            'step': {'__value': [1, 2]},
            'time': {'__value': [0.5, 1.0], '@unit': 'ps'},
            'value': {'__value': [3.0, 4.0], '@unit': 'kilojoule'},
        }

        result = parser.get_step_data(data, 2)

        assert result['time'].to('ps').magnitude == approx(1.0)
        assert result['value'].to('kilojoule').magnitude == approx(4.0)

    def test_maps_cell_data_for_a_trajectory_step(self):
        parser = H5MDH5Parser(
            data={
                'particles': {
                    'all': {
                        'box': {
                            '@boundary': [True, True, True],
                            'edges': {
                                'step': {'__value': [1, 2]},
                                'time': {'__value': [0.5, 1.0], '@unit': 'ps'},
                                'value': {
                                    '__value': np.array(
                                        [np.eye(3), 2 * np.eye(3)]
                                    ),
                                    '@unit': 'angstrom',
                                },
                            },
                        },
                    }
                }
            }
        )

        result = parser.get_cell_data({'step': 2})

        assert result['boundary'] == [True, True, True]
        assert_approx(result['lattice_vectors'].to('angstrom').magnitude, 2 * np.eye(3))

    def test_maps_species_labels_only_on_topology_frame(self):
        parser = H5MDH5Parser(
            data={'particles': {'all': {'species_label': ['H', 'O']}}}
        )

        result = parser.to_species_labels(
            {'step': 0, 'frame_index': 0}, path='particles.all.species_label'
        )

        assert [payload['chemical_symbol'] for payload in result] == ['H', 'O']
        assert parser.to_species_labels(
            {'step': 0, 'frame_index': 1}, path='particles.all.species_label'
        ) == []

    def test_maps_top_level_system_quantity(self):
        parser = H5MDH5Parser(data={'connectivity': {'bonds': [[0, 1]]}})

        assert parser.get_top_system_quantity(
            {'step': 0}, path='connectivity.bonds'
        ) == [[0, 1]]
        assert parser.get_top_system_quantity(
            {}, path='connectivity.bonds'
        ) == []

    def test_maps_trajectory_frames_and_filters_steps(self):
        parser = H5MDH5Parser()
        parser.trajectory_steps = [0, 2]
        source = {
            'position': {
                'step': {'__value': [0, 1, 2]},
                'time': {'__value': [0.0, 1.0, 2.0], '@unit': 'ps'},
                'value': {
                    '__value': np.arange(27).reshape(3, 3, 3),
                    '@unit': 'angstrom',
                },
            },
            'velocity': {
                'value': {
                    '__value': np.ones((3, 3, 3)),
                    '@unit': 'angstrom / ps',
                }
            },
        }

        frames = parser.get_traj_data(source)

        assert [frame['step'] for frame in frames] == [0, 2]
        assert (
            frames[1]['positions'].to('angstrom').magnitude[0, 0]
            == approx(18)
        )
        assert frames[0]['velocities'].to('angstrom / ps').magnitude.shape == (3, 3)

    def test_rejects_mismatched_trajectory_lengths(self):
        parser = H5MDH5Parser()
        source = {
            'position': {
                'step': {'__value': [0, 1]},
                'time': {'__value': [0.0]},
                'value': {'__value': np.zeros((2, 1, 3))},
            }
        }

        assert parser.get_traj_data(source) == []

    def test_collects_configurational_output_steps(self):
        parser = H5MDH5Parser()
        source = {
            'energies': {
                'total': {
                    '@type': 'configurational',
                    'step': {'__value': [1, 2]},
                    'time': {'__value': [0.5, 1.0], '@unit': 'ps'},
                }
            },
            'temperatures': {
                '@type': 'configurational',
                'step': {'__value': [2, 3]},
                'time': {'__value': [1.0, 1.5], '@unit': 'ps'},
            },
        }

        result = parser.get_output_steps(source)

        assert [item['step'] for item in result] == [1, 2, 3]
        assert [item['time'].to('ps').magnitude for item in result] == [0.5, 1.0, 1.5]

    def test_maps_observable_contributions_at_a_step(self):
        parser = H5MDH5Parser(
            data={
                'observables': {
                    'total_energy': {
                        'potential': {
                            'step': {'__value': [0, 2]},
                            'time': {'__value': [0.0, 2.0]},
                            'value': {
                                '__value': [1.0, 3.0],
                                '@unit': 'kilojoule',
                            },
                        }
                    }
                }
            }
        )

        result = parser.get_contributions(
            {'step': 2}, path='observables.total_energy'
        )

        assert result[0]['name'] == 'potential'
        assert result[0]['value'].to('kilojoule').magnitude == approx(3.0)

    def test_maps_ensemble_output_labels_and_units(self):
        parser = H5MDH5Parser()
        source = {
            'water': {
                'value': {'__value': np.array([1.0, 2.0]), '@unit': 'nm'},
                'n_times': {'__value': 2},
            }
        }

        result = parser.get_output_data(source)

        assert result[0]['label'] == 'water'
        assert_approx(result[0]['value'].to('nm').magnitude, [1.0, 2.0])

    def test_maps_configurational_output_at_a_step(self):
        parser = H5MDH5Parser(
            data={
                'observables': {
                    'temperatures': {
                        '@type': 'configurational',
                        'step': {'__value': [1, 2]},
                        'time': {'__value': [0.5, 1.0]},
                        'value': {'__value': [300.0, 310.0], '@unit': 'K'},
                    }
                }
            }
        )

        result = parser.get_configurational_output(
            {'step': 2}, path='observables.temperatures'
        )

        assert result.to('K').magnitude == approx(310.0)

    def test_maps_custom_configurational_outputs(self):
        parser = H5MDH5Parser(
            data={
                'observables': {
                    'custom_pressure': {
                        'step': {'__value': [1, 2]},
                        'time': {'__value': [0.5, 1.0]},
                        'value': {
                            '__value': [1.0, 2.0],
                            '@unit': 'bar',
                        },
                    },
                    'ignored': {'step': {'__value': [1]}, 'value': {'__value': [0]}},
                }
            }
        )

        result = parser.get_custom_outputs(
            {'step': 2}, path='observables', include=['custom_pressure']
        )

        assert result == [
            {
                'name': 'custom_pressure',
                'value': approx(2.0),
                'time': approx(1.0),
                'unit': 'bar',
            }
        ]

    def test_maps_custom_ensemble_outputs_by_type_and_exclusion(self):
        parser = H5MDH5Parser()
        source = {
            'rdf': {
                '@type': 'ensemble_average',
                'bins': {'__value': [0.1, 0.2], '@unit': 'nm'},
                'value': {'__value': [1.0, 2.0]},
            },
            'vacf': {
                '@type': 'correlation_function',
                'times': {'__value': [0.0, 1.0], '@unit': 'ps'},
                'value': {'__value': [1.0, 0.5], '@unit': 'nm ** 2 / ps ** 2'},
            },
            'standard': {'@type': 'ensemble_average', 'value': {'__value': [9.0]}},
        }

        result = parser.get_custom_ensemble_outputs(
            source, exclude=['standard'], type_filter=['ensemble_average']
        )

        assert len(result) == 1
        assert result[0]['label'] == 'rdf'
        assert_approx(result[0]['value_magnitude'], [1.0, 2.0])
        assert 'value_unit' not in result[0]

    def test_maps_direct_ensemble_output(self):
        parser = H5MDH5Parser()
        source = {
            'rdf': {
                'value': {'__value': [1.0, 2.0], '@unit': '1 / nm'},
                'bins': {'__value': [0.1, 0.2], '@unit': 'nm'},
            }
        }

        result = parser.get_ensemble_output(source)

        assert result[0]['label'] == 'rdf'
        assert_approx(result[0]['value'].to('1 / nm').magnitude, [1.0, 2.0])
