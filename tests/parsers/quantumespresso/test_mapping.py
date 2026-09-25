import numpy as np
import pytest
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Program, Simulation

from nomad_simulation_parsers.parsers.quantumespresso.gipaw.parser import (
    GIPAWMainfileTextParser,
    GIPAWMainfileXMLParser,
)
from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    MainfileTextParser,
    MainfileXMLParser,
    QuantumEspressoMetainfoParser,
    get_program_name_version,
)
from nomad_simulation_parsers.parsers.quantumespresso.phonon.parser import (
    PhononMainfileParser,
)
from nomad_simulation_parsers.parsers.quantumespresso.pwscf.parser import (
    PWSCFMainfileTextParser,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import common
from tests.parsers.common import approx, assert_approx


@pytest.fixture
def text_parser():
    return MainfileTextParser()


@pytest.fixture
def xml_parser():
    return MainfileXMLParser()


@pytest.fixture
def gipaw_text_parser():
    return GIPAWMainfileTextParser()


@pytest.fixture
def gipaw_xml_parser():
    return GIPAWMainfileXMLParser()


@pytest.fixture
def pwscf_text_parser():
    return PWSCFMainfileTextParser()


@pytest.mark.unit
class TestQuantumEspressoTextMapping:
    def test_common_and_module_mapping_keys_are_distinct(self):
        keys = {
            common.OUT_KEY,
            common.GIPAW_PROPERTIES_OUT_KEY,
            common.PWSCF_OUT_KEY,
        }

        assert len(keys) == 3

    def test_maps_text_helpers_and_topology_values(self, text_parser):
        parser = text_parser
        parser._data = {
            'header': {
                'alat': 2.0,
                'starting_magnetization': [0.1],
            }
        }

        assert parser.get_version(['PWSCF', 'v.7.3']) == '7.3'
        assert parser.get_datetime('21 Feb 2024 16:33:02').year == 2024
        assert parser.get_header('alat') == 2.0
        assert parser.get_n_spin_channels() == 2
        assert_approx(
            parser.get_value(
                {'positions': [[1.0, 0.0, 0.0]], 'units': 'alat'},
                key='positions',
            ),
            [[2.0, 0.0, 0.0]],
        )
        assert parser.get_topology_value({'frame_index': 1, 'x': [1]}) is None
        assert parser.get_periodic_boundary_conditions(
            {'simulation_cell': np.eye(3)}
        ) == [True, True, True]

    def test_maps_phonon_energy_with_common_outputs(self):
        source = PhononMainfileParser()
        source._data = {
            'header': {
                'program_name_version': ['PHONON', 'v.7.0'],
                'start_date_time': '1Jan2024 00:00:00',
            },
            'calculation': [
                {
                    'energies': {'energy_total': -10 * ureg.rydberg},
                    'simulation_cell': np.eye(3) * ureg.bohr,
                    'labels_positions': {
                        'labels': ['Si'],
                        'positions': np.zeros((1, 3)) * ureg.bohr,
                    },
                }
            ],
        }
        target = QuantumEspressoMetainfoParser(
            data_object=Simulation(program=Program(name='Quantum Espresso'))
        )
        target.annotation_key = common.OUT_KEY

        source.convert(target)

        outputs = target.data_object.outputs
        assert len(outputs) == 1
        assert outputs[0].m_def.qualified_name() == (
            'nomad_simulations.schema_packages.outputs.Outputs'
        )
        assert outputs[0].total_energies[0].value.to('rydberg').magnitude == approx(
            -10
        )

    def test_maps_program_name_and_version(self):
        assert get_program_name_version('Program PWSCF v.7.3 starts') == (
            'pwscf',
            (7, 3),
        )

    @pytest.mark.parametrize(
        ('source', 'expected'),
        [
            ('PBE', 'PBE'),
            ('SLA PW PBX PBC ( 1 4 3 4 0 0)', 'PBE'),
            ('SLA PW PSX PSC', 'PBEsol'),
        ],
    )
    def test_maps_xc_functional_keys(self, source, expected, text_parser):
        parser = text_parser
        assert parser.get_functional_key(source) == expected

    def test_maps_energy_contributions(self, text_parser):
        parser = text_parser
        source = {
            'energy_total': -10.0 * ureg.rydberg,
            'energy_total_one_electron': -8.0 * ureg.rydberg,
            'energy_total_hartree': -2.0 * ureg.rydberg,
        }

        contributions = parser.get_energy_contributions(source)

        assert [item['name'] for item in contributions] == [
            'one_electron',
            'hartree',
        ]
        assert contributions[0]['value'] == approx(-8.0)


@pytest.mark.unit
class TestQuantumEspressoXMLMapping:
    def test_maps_xml_helpers_and_forces(self, xml_parser):
        parser = xml_parser
        parser._data = {'@Units': 'Hartree atomic units'}

        assert parser.get_datetime('21Feb2024', '16:33:02').year == 2024
        assert parser.apply_unit(1.0, name='energy').to('hartree').magnitude == approx(
            1
        )
        assert_approx(parser.get_forces(np.arange(6.0)), [[0, 1, 2], [3, 4, 5]])
        assert parser.get_periodic_boundary_conditions(np.eye(3)) == [True] * 3

    def test_maps_energy_contributions_and_periodicity(self, xml_parser):
        parser = xml_parser
        source = {
            'etot': -1.0 * ureg.hartree,
            'one_electron': -0.8 * ureg.hartree,
        }

        contributions = parser.get_energy_contributions(source)

        assert contributions == [{'value': -0.8 * ureg.hartree, 'name': 'one_electron'}]
        assert parser.get_periodic_boundary_conditions(np.eye(3)) == [True] * 3


@pytest.mark.unit
class TestGIPAWMapping:
    def test_maps_text_magnetic_properties(self, gipaw_text_parser):
        parser = gipaw_text_parser
        source = {
            'ms_list': [
                ['Si', 1, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ],
            'chi_bare_pGv': np.diag([1.0, 2.0, 3.0]),
            'chi_bare_vGv': np.diag([3.0, 4.0, 5.0]),
        }

        mapped = parser.get_gipaw_text(source)[0]

        assert_approx(
            mapped['magnetic_shieldings'][0]['value'].magnitude,
            np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) * 1e-6,
        )
        assert_approx(
            mapped['magnetic_susceptibilities']['value'],
            np.diag([2.0, 3.0, 4.0]),
        )
        assert_approx(
            mapped['magnetic_susceptibilities']['value_pgv_approx'],
            np.diag([1.0, 2.0, 3.0]),
        )

    @pytest.mark.parametrize(
        ('job', 'method', 'expected'),
        [
            ('nmr', 'get_magnetic_shieldings', np.eye(3) * 1e-6),
            ('efg', 'get_efg', np.eye(3)),
            ('hyperfine', 'get_hyperfine_dipolar', np.eye(3)),
            ('g-tensor', 'get_delta_g', np.eye(3)),
        ],
    )
    def test_maps_job_specific_tensors(self, job, method, expected, gipaw_xml_parser):
        parser = gipaw_xml_parser
        parser._data = {'input': {'job': job}}
        value = getattr(parser, method)({'__value': np.eye(3).reshape(-1)})

        assert_approx(
            value.magnitude if hasattr(value, 'magnitude') else value,
            expected,
        )

    def test_maps_hyperfine_fermi_contact_and_job_guards(self, gipaw_xml_parser):
        parser = gipaw_xml_parser
        parser._data = {'input': {'job': 'hyperfine'}}

        assert parser.get_hyperfine_fermi_contact({'__value': 2.5}) == approx(2.5)
        assert parser.get_magnetic_shieldings({'__value': np.eye(3)}) is None


@pytest.mark.unit
class TestPWSCFMapping:
    def test_maps_forces_and_configurations(self, pwscf_text_parser):
        parser = pwscf_text_parser
        force = np.arange(3.0)
        source = {
            'self_consistent': {
                'self_consistent': [
                    {'forces': force},
                    {'forces': force + 1},
                ]
            },
            'bfgs_geometry_optimization': {'self_consistent': [{'forces': force + 2}]},
        }

        configurations = parser.get_configurations(source)
        configuration_forces = parser.get_configuration_forces(source)

        assert len(configurations) == 3
        assert [config['frame_index'] for config in configurations] == [0, 1, 2]
        assert len(configuration_forces) == 2
        assert len(configuration_forces[0]) == 2
        assert_approx(configuration_forces[1][0], force + 2)
        contributions = parser.get_force_contributions({'forces_dispersion': force})
        assert contributions[0]['name'] == 'dispersion'

    def test_maps_eigenvalues_and_occupations(self, pwscf_text_parser):
        parser = pwscf_text_parser
        parser._data = {'header': {}}
        source = {
            'band_energies': [[-1.0, 0.5, 2.0]],
            'occupation_numbers': [[2.0, 1.0, 0.0]],
        }

        eigenvalues = parser.get_eigenvalues(source)

        assert len(eigenvalues) == 1
        assert eigenvalues[0]['n_levels'] == 3
        assert_approx(
            eigenvalues[0]['eigenvalues'].to('eV').magnitude,
            [[-1, 0.5, 2]],
        )
        assert_approx(eigenvalues[0]['occupations'], [[2, 1, 0]])

    def test_maps_scf_steps_and_reference_energy(self, pwscf_text_parser):
        parser = pwscf_text_parser
        source = {
            'iteration': [
                {
                    'energies': {
                        'energy_total': -10 * ureg.rydberg,
                        'energy_total_accuracy_estimate': 0.1 * ureg.rydberg,
                    },
                    'time': 1.5,
                    'threshold': 1e-4,
                    'ddv_scf': 2e-5,
                },
                {
                    'energies': {
                        'energy_total': -10.1 * ureg.rydberg,
                        'energy_total_accuracy_estimate': 0.01 * ureg.rydberg,
                    },
                    'time': 1.0,
                    'threshold': 1e-6,
                    'ddv_scf': 3e-7,
                },
            ],
            'homo_lumo': [-2.5],
        }

        steps = parser.get_scf_steps(source)

        assert_approx(
            [value.to('rydberg').magnitude for value in steps['energies_total']],
            [-10, -10.1],
        )
        assert_approx(
            [value.to('rydberg').magnitude for value in steps['delta_energies_total']],
            [0.1, 0.01],
        )
        assert steps['durations'] == [1.5, 1.0]
        assert parser.get_reference_energy(source).to('eV').magnitude == approx(-2.5)

    def test_maps_band_structures_and_single_point_workflow(self, pwscf_text_parser):
        parser = pwscf_text_parser
        parser._data = {'header': {'scf_threshold_energy_change': 1e-8 * ureg.rydberg}}
        source = {
            'band_energies': [[-1.0, 0.5, 2.0]],
            'occupation_numbers': [[2.0, 1.0, 0.0]],
            'homo_lumo': [-2.5],
        }

        bands = parser.get_band_structures(source)
        workflow = parser.build_workflow()

        assert len(bands) == 1
        assert_approx(bands[0]['value'].to('eV').magnitude, [[-1, 0.5, 2]])
        assert bands[0]['highest_occupied'].to('eV').magnitude == approx(-2.5)
        assert workflow.m_def.name == 'SinglePoint'
        assert len(workflow.method.convergence_targets) == 1
