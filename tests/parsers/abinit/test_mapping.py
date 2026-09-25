import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.outputs import Outputs

from nomad_simulation_parsers.parsers.abinit.parser import (
    AbinitArchiveWriter,
    AbinitMetainfoParser,
    DosParser,
    MainfileParser,
)
from nomad_simulation_parsers.schema_packages import abinit


class ParsedSource(dict):
    """Controlled output of AbinitOutParser used at the mapping boundary."""

    def __init__(self, values, input_vars):
        super().__init__(values)
        self.input_vars = input_vars


def source_parser(source):
    parser = MainfileParser()
    # MappingParser consumes ``data`` while MainfileParser helper methods consume
    # ``data_object``. Point both boundaries at the same controlled source.
    parser._data = source
    parser._data_object = source
    return parser


def controlled_source(**input_overrides):
    input_vars = {
        'natom': [[2]],
        'xcart': [np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]) * ureg.bohr],
        'znucl': [[14]],
        'typat': [[1, 1]],
        'ixc': [[11]],
        'nsppol': [[1]],
    }
    input_vars.update(input_overrides)
    values = {
        'program_version': '9.10.4',
        'x_abinit_start_date': 'Mon 01 Jan 2024',
        'x_abinit_start_time': '12h30',
        'dataset': [
            {
                'x_abinit_vprim': np.eye(3) * 4.0,
                'results': {
                    'energy_total': -3.5 * ureg.hartree,
                    'cartesian_forces': np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
                    * ureg.hartree
                    / ureg.bohr,
                    # Rows contain k-point metadata followed by band energies.
                    'eigenvalues': np.array(
                        [
                            [1.0, 2.0, 0.5, 0.0, 0.0, 0.0, -1.0, 1.0],
                            [2.0, 2.0, 0.5, 0.5, 0.0, 0.0, -0.5, 2.0],
                        ]
                    ),
                    'occupation_numbers': np.array([[2.0, 0.0], [2.0, 0.0]]),
                },
                'self_consistent': {
                    'energy_total_scf_iteration': [
                        [-3.0, 0.1],
                        [-3.5, 0.0],
                    ]
                },
            }
        ],
    }
    return ParsedSource(values, input_vars)


@pytest.mark.unit
class TestOutMapping:
    def test_maps_controlled_source_to_metainfo_with_units_and_ordering(self):
        parser = source_parser(controlled_source())
        archive = EntryArchive()
        archive.data = Simulation(program=Program(name='ABINIT'))
        target = AbinitMetainfoParser()
        target.annotation_key = abinit.OUT_KEY
        target.data_object = archive.data

        parser.convert(target)

        simulation = archive.data
        assert simulation.program.name == 'ABINIT'
        assert simulation.program.version == '9.10.4'
        assert simulation.datetime.isoformat() == '2024-01-01T12:30:00+00:00'

        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Si',
            'Si',
        ]
        np.testing.assert_allclose(
            system.positions.to('bohr').magnitude,
            [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        )
        np.testing.assert_allclose(
            system.lattice_vectors.to('bohr').magnitude, np.eye(3) * 4.0
        )
        assert list(system.periodic_boundary_conditions) == [True, True, True]

        components = simulation.model_method[0].xc.components
        assert [component.canonical_label for component in components] == [
            'GGA_X_PBE',
            'GGA_C_PBE',
        ]

        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert output.total_energies[0].value.to('hartree').magnitude == pytest.approx(
            -3.5
        )
        np.testing.assert_allclose(
            output.total_forces[0].value.to('hartree / bohr').magnitude,
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        )
        np.testing.assert_allclose(
            output.scf_steps.energies_total.to('hartree').magnitude, [-3.0, -3.5]
        )
        # Zero is data and must not be confused with an absent delta-energy field.
        np.testing.assert_allclose(
            output.scf_steps.delta_energies_total.to('hartree').magnitude, [0.1, 0.0]
        )

        assert len(output.electronic_band_structures) == 1
        band_structure = output.electronic_band_structures[0]
        assert band_structure.value.shape == (2, 2)
        np.testing.assert_allclose(
            band_structure.value.to('hartree').magnitude,
            [[-1.0, 1.0], [-0.5, 2.0]],
        )
        np.testing.assert_allclose(band_structure.occupation, [[2.0, 0.0], [2.0, 0.0]])
        assert band_structure.spin_channel is None
        assert output.electronic_band_gaps[0].value.to(
            'hartree'
        ).magnitude == pytest.approx(1.5)

    def test_missing_scf_iterations_map_to_no_scf_section(self):
        source = controlled_source()
        source['dataset'][0].pop('self_consistent')
        parser = source_parser(source)

        assert parser.get_scf_steps(source['dataset'][0]) == {}
        assert 'scf_steps' not in parser.get_outputs()[0]


@pytest.mark.unit
class TestDosMapping:
    def test_maps_spin_resolved_dos_with_energy_and_density_units(self):
        parser = DosParser()
        parser._data = {
            'nspinpol': 2,
            'data': np.array([[-1.0, 10.0], [0.0, 20.0], [-1.0, 30.0], [0.0, 40.0]]),
        }
        archive = EntryArchive()
        archive.data = Simulation(outputs=[Outputs()])
        target = AbinitMetainfoParser()
        target.annotation_key = abinit.DOS_KEY
        target.data_object = archive.data

        parser.convert(target, update_mode='merge@-1')

        dos = archive.data.outputs[0].electronic_dos
        assert len(dos) == 2
        np.testing.assert_allclose(
            dos[0].energies.points.to('hartree').magnitude, [-1.0, 0.0]
        )
        np.testing.assert_allclose(
            dos[0].value.to('1 / hartree').magnitude, [10.0, 20.0]
        )
        np.testing.assert_allclose(
            dos[1].value.to('1 / hartree').magnitude, [30.0, 40.0]
        )


@pytest.mark.unit
class TestWorkflowMapping:
    def test_maps_geometry_workflow_and_convergence_contract(self):
        source = controlled_source(
            ionmov=[[2]],
            vis=[[100.0]],
            optcell=[[0]],
            tolmxde=[[1.0e-6]],
            tolmxf=[[2.0e-4]],
        )
        writer = AbinitArchiveWriter()
        writer.archive = EntryArchive()
        writer.mainfile_parser = source_parser(source)

        writer.parse_workflow()

        workflow = writer.archive.workflow2
        assert workflow.m_def.name == 'GeometryOptimization'
        assert workflow.method.optimization_type == 'atomic'
        assert workflow.method.optimization_method == 'bfgs'
        targets = {
            target.m_def.name: target for target in workflow.method.convergence_targets
        }
        assert targets['EnergyConvergenceTarget'].threshold.to(
            'hartree'
        ).magnitude == pytest.approx(1.0e-6)
        assert targets['ForceConvergenceTarget'].threshold.to(
            'hartree / bohr'
        ).magnitude == pytest.approx(2.0e-4)

    def test_maps_single_point_workflow_without_optional_convergence(self):
        source = controlled_source(ionmov=[[0]], vis=[[100.0]])
        writer = AbinitArchiveWriter()
        writer.archive = EntryArchive()
        writer.mainfile_parser = source_parser(source)

        writer.parse_workflow()

        assert writer.archive.workflow2.m_def.name == 'SinglePoint'
        assert writer.archive.workflow2.method is None
