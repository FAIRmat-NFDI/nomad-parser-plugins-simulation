import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.ams.parser import MainfileParser


@pytest.mark.unit
class TestOutMapping:
    def test_maps_periodic_boundary_conditions_from_lattice_shape(self):
        parser = MainfileParser()

        assert parser.get_periodic_boundary_conditions({'lattice_vectors': None}) == [
            False,
            False,
            False,
        ]
        assert parser.get_periodic_boundary_conditions(
            np.eye(2, 3) * ureg.angstrom
        ) == [
            True,
            True,
            False,
        ]
        assert parser.get_periodic_boundary_conditions(np.eye(3) * ureg.angstrom) == [
            True,
            True,
            True,
        ]


@pytest.mark.unit
class TestWorkflowMapping:
    def test_builds_workflows_and_convergence_targets(self):
        parser = MainfileParser()

        workflow = parser.build_workflow(
            {
                'geometry_optimization': {
                    'convergence_tolerance_force_maximum': 2e-4
                    * ureg.hartree
                    / ureg.bohr,
                    'convergence_tolerance_energy_difference': 1e-6 * ureg.hartree,
                }
            }
        )

        assert workflow.m_def.name == 'GeometryOptimization'
        assert {
            target.m_def.name for target in workflow.method.convergence_targets
        } == {
            'EnergyConvergenceTarget',
            'ForceConvergenceTarget',
        }
        assert parser.build_workflow({'molecular_dynamics': {}}).m_def.name == (
            'MolecularDynamics'
        )


@pytest.mark.unit
class TestDosMapping:
    def test_extracts_dos_columns(self):
        parser = MainfileParser()
        result = parser.get_dos({'dos': np.array([[-1.0, 2.0, 3.0], [0.0, 4.0, 5.0]])})

        assert len(result) == 2
        np.testing.assert_allclose(result[0]['energies'], [-1.0, 0.0])
        np.testing.assert_allclose(result[0]['value'], [2.0, 4.0])
        np.testing.assert_allclose(result[1]['value'], [3.0, 5.0])
