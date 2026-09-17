from types import SimpleNamespace

import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.gpaw.parser import GPWParser
from tests.parsers.common import approx, assert_approx


class ControlledSource:
    def __init__(self):
        self.parameters = {
            'energy_total': -3.5,
            'energy_XC': -1.0,
            'energyerror': 1.0e-6,
            'fermilevel': 0.0,
            'converged': True,
        }
        self.arrays = {
            'atom_forces_free': np.array([[1.0, 2.0, 3.0]]),
            'eigenvalues': np.array([[-1.0, 0.5]]),
            'occupation': np.array([[1.0, 0.0]]),
        }

    def get_parameter(self, key):
        return self.parameters.get(key)

    def get_array(self, key):
        return self.arrays.get(key)


def mapper(source=None):
    source = source or ControlledSource()

    def apply_unit(value, unit):
        if value is None:
            return None
        unit_name = {'energyunit': 'eV', 'lengthunit': 'angstrom'}.get(unit, unit)
        return value * getattr(ureg, unit_name)

    data_object = SimpleNamespace(
        parser=source,
        apply_unit=apply_unit,
    )
    parser = GPWParser()
    parser.data_object = data_object
    return parser


@pytest.mark.unit
class TestGPAWMapping:
    def test_maps_energies_forces_and_scf_state(self):
        parser = mapper()

        energies = parser.get_energies()
        assert energies['total'].to('eV').magnitude == approx(-3.5)
        assert energies['contributions'][0]['name'] == 'XC'

        forces = parser.get_forces()
        assert_approx(forces['value'].to('eV / angstrom').magnitude, [[1.0, 2.0, 3.0]])
        assert parser.get_scf_steps()['code_specific_quantities'] == {
            'converged': True,
            'energyerror': 1.0e-6,
        }

    def test_maps_eigenvalues_occupations_and_gap(self):
        parser = mapper()

        eigenvalues = parser.get_eigenvalues()[0]
        assert_approx(eigenvalues['value'].to('eV').magnitude, [[-1.0, 0.5]])
        assert_approx(eigenvalues['occupation'], [[1.0, 0.0]])
        assert eigenvalues['n_levels'] == 2
        assert eigenvalues['highest_occupied'].to('eV').magnitude == approx(0.0)

        gap = parser.get_band_gaps()[0]
        assert gap['value'].to('eV').magnitude == approx(1.5)

    def test_builds_single_point_workflow_with_energy_target(self):
        workflow = mapper().build_workflow()

        assert workflow.m_def.name == 'SinglePoint'
        target = workflow.method.convergence_targets[0]
        assert target.m_def.name == 'EnergyConvergenceTarget'
        assert target.threshold.to('eV').magnitude == approx(1.0e-6)
