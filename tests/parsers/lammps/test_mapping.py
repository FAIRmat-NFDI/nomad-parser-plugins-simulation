"""Isolated LAMMPS method-to-schema mapping contracts."""

import numpy as np
import pytest

from nomad_simulation_parsers.parsers.lammps.file_parsers import get_unit
from nomad_simulation_parsers.parsers.lammps.md_workflow_method import build_md_method
from tests.parsers.common import approx


@pytest.mark.unit
def test_maps_controlled_npt_settings_without_reading_a_log():
    settings = {
        'fixes': [
            [
                '1',
                'all',
                'npt',
                'temp',
                '300',
                '310',
                '100',
                'iso',
                '1',
                '2',
                '1000',
            ]
        ],
        'dumps': [
            [
                '1',
                'all',
                'custom',
                '50',
                'trajectory.dump',
                'id',
                'x',
                'y',
                'z',
                'vx',
                'vy',
                'vz',
                50,
            ]
        ],
        'timestep': 2.0,
        'thermo': 100,
        'run_style': 'verlet',
        'n_steps': 500,
    }

    method = build_md_method(settings, get_unit('real'))

    assert method.m_def.name == 'MolecularDynamicsMethod'
    assert method.thermodynamic_ensemble == 'NPT'
    assert method.integrator_type == 'velocity_verlet'
    assert method.integration_timestep.to('fs').magnitude == approx(2)
    assert method.n_steps == 500
    assert method.coordinate_save_frequency == 50
    assert method.velocity_save_frequency == 50
    assert method.force_save_frequency is None
    assert method.thermodynamics_save_frequency == 100
    thermostat = method.thermostat_parameters[0]
    assert thermostat.temperature_profile == 'linear'
    assert thermostat.reference_temperature_start.to('K').magnitude == approx(300)
    assert thermostat.reference_temperature_end.to('K').magnitude == approx(310)
    assert thermostat.coupling_constant.to('fs').magnitude == approx(100)
    barostat = method.barostat_parameters[0]
    assert barostat.coupling_type == 'isotropic'
    assert barostat.pressure_profile == 'linear'
    np.testing.assert_allclose(
        barostat.reference_pressure_start.to('atm').magnitude,
        np.eye(3),
    )
    np.testing.assert_allclose(
        barostat.reference_pressure_end.to('atm').magnitude,
        np.eye(3) * 2,
    )
    np.testing.assert_allclose(
        barostat.coupling_constant.to('fs').magnitude,
        np.eye(3) * 1000,
    )


@pytest.mark.unit
def test_zero_steps_are_preserved_without_inventing_missing_settings():
    settings = {
        'fixes': [],
        'dumps': [],
        'timestep': None,
        'thermo': None,
        'run_style': None,
        'n_steps': 0,
    }

    method = build_md_method(settings, {})

    assert method.n_steps == 0
    assert method.m_to_dict()['n_steps'] == 0
    assert method.thermodynamic_ensemble is None
    assert method.integration_timestep is None
    assert method.thermodynamics_save_frequency is None
    assert len(method.thermostat_parameters) == 0
    assert len(method.barostat_parameters) == 0
