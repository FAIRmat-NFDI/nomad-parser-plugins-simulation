from types import SimpleNamespace

import numpy as np
import phonopy
import pytest
from nomad_file_parser.mapping_parser import MetainfoParser

from nomad_simulation_parsers.parsers.phonopy import parser as parser_module
from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyYamlParser
from nomad_simulation_parsers.schema_packages.phonopy import (
    PHONOPY_KEY,
    PhonopySimulation,
)


@pytest.fixture
def source_parser():
    unitcell = SimpleNamespace(
        cell=np.eye(3) * 2,
        positions=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
        symbols=['Si', 'Si'],
    )
    supercell = SimpleNamespace(
        cell=np.eye(3) * 4,
        positions=np.zeros((2, 3)),
        symbols=['Si', 'Si'],
    )
    phonopy_obj = SimpleNamespace(
        unitcell=unitcell,
        supercell=supercell,
        supercell_matrix=np.diag([2, 2, 2]),
        displacements=np.array([[[0.0, 0.0, 0.0], [0.01, 0.0, 0.0]]]),
        symmetry=SimpleNamespace(tolerance=1e-5),
        calculator='vasp',
        nac_params=None,
        force_constants=np.ones((2, 2, 3, 3)),
    )
    parser = PhonopyYamlParser(data_object=phonopy_obj)
    parser._phonon_properties = SimpleNamespace(
        mesh=[2, 2, 2], frequencies=np.array([-1.0, 0.0, 1.0])
    )
    return parser


@pytest.mark.unit
class TestPhonopyYamlMapping:
    def test_maps_program_systems_and_method(self, source_parser):
        assert source_parser.get_program() == {
            'name': 'Phonopy',
            'version': phonopy.__version__,
        }

        systems = source_parser.get_model_system()
        assert len(systems) == 2
        assert [
            state['chemical_symbol'] for state in systems[0]['particle_states']
        ] == [
            'Si',
            'Si',
        ]
        np.testing.assert_allclose(
            systems[0]['positions'].to('angstrom').magnitude,
            [[0, 0, 0], [1, 1, 1]],
        )
        np.testing.assert_allclose(
            systems[1]['representations'][0]['lattice_vectors']
            .to('angstrom')
            .magnitude,
            np.eye(3) * 4,
        )
        np.testing.assert_array_equal(
            systems[1]['representations'][0]['supercell_matrix'],
            np.diag([2, 2, 2]),
        )
        assert systems[0]['representations'][0]['periodic_boundary_conditions'] == [
            True,
            True,
            True,
        ]

        method = source_parser.get_model_method()[0]
        assert method['name'] == 'harmonic lattice dynamics'
        assert method['type'] == 'finite displacement'
        assert method['force_calculator'] == 'vasp'
        assert method['with_non_analytic_correction'] is False
        assert method['displacement'].to('angstrom').magnitude == pytest.approx(0.01)
        assert method['symmetry_tolerance'].to('angstrom').magnitude == pytest.approx(
            1e-5
        )
        assert method['mesh_density'].to(
            '1 / angstrom ** 3'
        ).magnitude == pytest.approx(1.0)

    def test_maps_force_constants_and_calculated_outputs(
        self, source_parser, monkeypatch
    ):
        monkeypatch.setattr(
            parser_module,
            'get_bandstructures',
            lambda _: [
                {
                    'frequencies': np.array([[1.0, 2.0]]),
                    'kpoints': np.array([[0.0, 0.0, 0.0]]),
                    'labels': ['Γ', 'X'],
                }
            ],
        )
        monkeypatch.setattr(
            parser_module,
            'get_dos',
            lambda _: [
                {'frequencies': np.array([1.0, 2.0]), 'dos': np.array([3.0, 4.0])}
            ],
        )
        monkeypatch.setattr(
            parser_module,
            'get_thermodynamic_properties',
            lambda _: [
                {'temperature': 100.0, 'free_energy': 2.0, 'heat_capacity': 3.0},
                {'temperature': 300.0, 'free_energy': 4.0, 'heat_capacity': 5.0},
            ],
        )

        force_constants = source_parser.get_force_constants()[0]['value']
        np.testing.assert_allclose(
            force_constants.to('eV / angstrom ** 2').magnitude,
            np.ones((2, 2, 3, 3)),
        )
        bands = source_parser.get_phonon_band_structures()[0]
        assert bands['n_bands'] == 2
        assert bands['endpoints_labels'] == ['Γ', 'X']
        np.testing.assert_allclose(bands['value'].to('joule').magnitude, [[1, 2]])
        np.testing.assert_allclose(bands['q_path']['points'], [[0, 0, 0]])
        dos = source_parser.get_phonon_dos()[0]
        np.testing.assert_allclose(dos['value'].to('1 / joule').magnitude, [3, 4])
        np.testing.assert_allclose(
            dos['frequencies']['points'].to('joule').magnitude, [1, 2]
        )
        thermal = source_parser.get_vibrational_thermodynamics()
        np.testing.assert_allclose(
            thermal['vibrational_free_energies'][0]['value'].to('joule').magnitude,
            [2, 4],
        )
        np.testing.assert_allclose(
            thermal['vibrational_heat_capacities'][0]['value']
            .to('joule / kelvin')
            .magnitude,
            [3, 5],
        )
        np.testing.assert_allclose(
            thermal['vibrational_free_energies'][0]['temperatures']['points']
            .to('kelvin')
            .magnitude,
            [100, 300],
        )
        assert source_parser.get_vibrational_thermodynamics() is thermal
        assert source_parser.get_imaginary_frequencies() == 1

    def test_annotations_convert_to_phonopy_simulation(
        self, source_parser, monkeypatch
    ):
        monkeypatch.setattr(parser_module, 'get_bandstructures', lambda _: [])
        monkeypatch.setattr(parser_module, 'get_dos', lambda _: [])
        monkeypatch.setattr(parser_module, 'get_thermodynamic_properties', lambda _: [])
        target = MetainfoParser(data_object=PhonopySimulation())
        target.annotation_key = PHONOPY_KEY

        source_parser.convert(target)

        simulation = target.data_object
        assert simulation.program.name == 'Phonopy'
        assert len(simulation.model_system) == 2
        assert len(simulation.model_method) == 1
        assert simulation.model_method[0].displacement.to('angstrom').magnitude == (
            pytest.approx(0.01)
        )
        assert len(simulation.outputs) == 1
        np.testing.assert_allclose(
            simulation.outputs[0]
            .force_constants[0]
            .value.to('eV / angstrom ** 2')
            .magnitude,
            np.ones((2, 2, 3, 3)),
        )
        assert simulation.outputs[0].n_imaginary_frequencies == 1
