"""LAMMPS workflow method regressions, including command order and units."""

import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.file_parsers import get_unit
from nomad_simulation_parsers.parsers.lammps.md_workflow_method import (
    build_md_method,
    selected_run,
)
from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser
from tests.parsers.common import approx


@pytest.mark.unit
def test_uses_settings_active_at_last_run():
    settings = selected_run(
        [
            'units real\n',
            'fix 1 all nvt temp 300 300 100\n',
            'dump 1 all custom 100 out.xyz id x y z vx vy vz\n',
            'run 100\n',
            'unfix 1\n',
            'fix 2 all nve\n',
            'reset_timestep 500\n',
            'dump_modify 1 every 20\n',
            'run 550 upto\n',
        ]
    )
    method = build_md_method(settings, get_unit('real'))

    assert method.thermodynamic_ensemble == 'NVE'
    assert method.n_steps == 50
    assert method.coordinate_save_frequency == 20
    assert method.velocity_save_frequency == 20
    assert len(method.thermostat_parameters) == 0


@pytest.mark.unit
def test_ramped_nph_pressure_and_respa():
    settings = selected_run(
        [
            'run_style respa 2 2 bond 1 pair 2\n',
            'fix baro all nph iso 1 2 1000\n',
            'timestep 0.5\n',
            'run 25\n',
        ]
    )
    method = build_md_method(settings, get_unit('metal'))

    assert method.thermodynamic_ensemble == 'NPH'
    assert method.integrator_type == 'rRESPA_multitimescale'
    assert method.integration_timestep.to('ps').magnitude == approx(0.5)
    barostat = method.barostat_parameters[0]
    assert barostat.pressure_profile == 'linear'
    assert np.diag(barostat.reference_pressure_start.to('bar').magnitude) == approx(
        [1, 1, 1]
    )
    assert np.diag(barostat.reference_pressure_end.to('bar').magnitude) == approx(
        [2, 2, 2]
    )
    assert np.diag(barostat.coupling_constant.to('ps').magnitude) == approx(
        [1000, 1000, 1000]
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ('control', 'expected_type'),
    [
        ('fix 2 all langevin 250 300 50 12345', 'langevin_schneider'),
        ('fix 2 all temp/berendsen 250 300 50', 'berendsen'),
    ],
)
def test_separate_temperature_controls(control, expected_type):
    settings = selected_run(['fix 1 all nve\n', f'{control}\n', 'run 100\n'])
    method = build_md_method(settings, get_unit('real'))

    assert method.thermodynamic_ensemble == 'NVT'
    thermostat = method.thermostat_parameters[0]
    assert thermostat.thermostat_type == expected_type
    assert thermostat.temperature_profile == 'linear'
    assert thermostat.reference_temperature_start.to('K').magnitude == approx(250)
    assert thermostat.reference_temperature_end.to('K').magnitude == approx(300)
    assert thermostat.coupling_constant.to('fs').magnitude == approx(50)


@pytest.mark.unit
@pytest.mark.parametrize(
    ('style', 'expected_type'),
    [
        ('temp/csvr', 'velocity_rescaling'),
        ('temp/csvr/kk', 'velocity_rescaling'),
        ('temp/csld', 'velocity_rescaling_langevin'),
    ],
)
def test_canonical_sampling_thermostats(style, expected_type):
    settings = selected_run(
        ['fix md all nve\n', f'fix bath all {style} 250 300 50 12345\n', 'run 100\n']
    )
    method = build_md_method(settings, get_unit('real'))

    assert method.thermodynamic_ensemble == 'NVT'
    thermostat = method.thermostat_parameters[0]
    assert thermostat.thermostat_type == expected_type
    assert thermostat.reference_temperature_start.to('K').magnitude == approx(250)
    assert thermostat.reference_temperature_end.to('K').magnitude == approx(300)
    assert thermostat.coupling_constant.to('fs').magnitude == approx(50)


@pytest.mark.unit
def test_berendsen_pressure_control_with_separate_integrator_and_thermostat():
    settings = selected_run(
        [
            'fix md all nve\n',
            'fix bath all temp/csvr 300 300 100 12345\n',
            'fix baro all press/berendsen iso 1 2 1000 modulus 10\n',
            'run 100\n',
        ]
    )
    method = build_md_method(settings, get_unit('metal'))

    assert method.thermodynamic_ensemble == 'NPT'
    assert method.integrator_type == 'velocity_verlet'
    barostat = method.barostat_parameters[0]
    assert barostat.barostat_type == 'berendsen'
    assert barostat.coupling_type == 'isotropic'
    assert barostat.pressure_profile == 'linear'
    assert np.diag(barostat.reference_pressure_start.to('bar').magnitude) == approx(
        [1, 1, 1]
    )
    assert np.diag(barostat.reference_pressure_end.to('bar').magnitude) == approx(
        [2, 2, 2]
    )
    assert np.diag(barostat.coupling_constant.to('ps').magnitude) == approx(
        [1000, 1000, 1000]
    )
    assert np.diag(barostat.compressibility.to('1/bar').magnitude) == approx(
        [0.1, 0.1, 0.1]
    )
    assert 'compressibility' in method.m_to_dict()['barostat_parameters'][0]


@pytest.mark.unit
def test_semi_isotropic_npt_pressure_targets():
    settings = selected_run(
        [
            'fix baro all npt temp 300 300 100 '
            'x 1 1 500 y 1 1 500 z 2 3 700 couple xy\n',
            'run 25\n',
        ]
    )
    method = build_md_method(settings, get_unit('metal'))

    assert method.thermodynamic_ensemble == 'NPT'
    barostat = method.barostat_parameters[0]
    assert barostat.coupling_type == 'semi_isotropic'
    assert barostat.pressure_profile == 'linear'
    assert np.diag(barostat.reference_pressure_start.to('bar').magnitude) == approx(
        [1, 1, 2]
    )
    assert np.diag(barostat.reference_pressure_end.to('bar').magnitude) == approx(
        [1, 1, 3]
    )
    assert np.diag(barostat.coupling_constant.to('ps').magnitude) == approx(
        [500, 500, 700]
    )


@pytest.mark.unit
def test_triclinic_nph_pressure_and_partial_targets():
    settings = selected_run(
        ['fix baro all nph tri 1 1 1000 xy 0.1 0.1 200\n', 'run 25\n']
    )
    method = build_md_method(settings, get_unit('metal'))

    barostat = method.barostat_parameters[0]
    assert barostat.coupling_type == 'anisotropic'
    assert barostat.reference_pressure.to('bar').magnitude[0, 1] == approx(0.1)
    assert barostat.reference_pressure.to('bar').magnitude[1, 0] == approx(0.1)
    assert barostat.coupling_constant.to('ps').magnitude[0, 1] == approx(200)

    partial = selected_run(['fix baro all nph x 1 1 500\n', 'run 25\n'])
    partial_method = build_md_method(partial, get_unit('metal'))
    partial_barostat = partial_method.barostat_parameters[0]
    assert partial_barostat.reference_pressure is None
    assert np.diag(partial_barostat.coupling_constant.to('ps').magnitude) == approx(
        [500, 0, 0]
    )


@pytest.mark.unit
def test_brownian_temperature_is_not_a_damping_time():
    settings = selected_run(['fix bd all brownian 300 12345 gamma_t 1\n', 'run 100\n'])
    method = build_md_method(settings, get_unit('real'))

    assert method.thermodynamic_ensemble == 'NVT'
    assert method.integrator_type == 'brownian'
    thermostat = method.thermostat_parameters[0]
    assert thermostat.thermostat_type == 'brownian'
    assert thermostat.reference_temperature.to('K').magnitude == approx(300)
    assert thermostat.coupling_constant is None


@pytest.mark.unit
def test_accelerated_nvt_style_keeps_the_same_method_semantics():
    settings = selected_run(['fix md all nvt/kk temp 300 300 100\n', 'run 10\n'])
    method = build_md_method(settings, get_unit('real'))

    assert method.thermodynamic_ensemble == 'NVT'
    assert method.integrator_type == 'velocity_verlet'
    assert method.thermostat_parameters[0].thermostat_type == 'nose_hoover'


@pytest.mark.unit
def test_group_specific_fixes_do_not_invent_a_global_ensemble():
    settings = selected_run(
        [
            'fix 1 mobile nvt temp 300 300 100\n',
            'fix 2 framework nve\n',
            'run 10\n',
        ]
    )
    method = build_md_method(settings, get_unit('real'))

    assert method.n_steps == 10
    assert method.thermodynamic_ensemble is None
    assert method.integrator_type is None
    assert len(method.thermostat_parameters) == 0


@pytest.mark.unit
def test_incomplete_nvt_control_does_not_crash_or_invent_values():
    settings = selected_run(['fix md all nvt\n', 'run 10\n'])
    method = build_md_method(settings, get_unit('real'))

    assert method.thermodynamic_ensemble == 'NVT'
    assert method.thermostat_parameters[0].thermostat_type == 'nose_hoover'
    assert method.thermostat_parameters[0].reference_temperature is None


@pytest.mark.unit
def test_ambiguous_dump_frequencies_and_lj_units_are_not_invented():
    settings = selected_run(
        [
            'fix 1 all nvt temp 8 8 0.1\n',
            'timestep 0.001\n',
            'dump 1 all xyz 10 a.xyz\n',
            'dump 2 all xyz 20 b.xyz\n',
            'run 1000\n',
        ]
    )
    method = build_md_method(settings, {})

    assert method.thermodynamic_ensemble == 'NVT'
    assert method.coordinate_save_frequency is None
    assert method.integration_timestep is None
    assert method.thermostat_parameters[0].reference_temperature is None


@pytest.mark.unit
def test_custom_dump_save_frequencies_are_distinct_from_thermo():
    settings = selected_run(
        [
            'fix md all nve\n',
            'dump d all custom 20 trajectory.dump id x y z vx vy vz fx fy fz\n',
            'thermo 100\n',
            'run 200\n',
        ]
    )
    method = build_md_method(settings, get_unit('real'))

    assert method.coordinate_save_frequency == 20
    assert method.velocity_save_frequency == 20
    assert method.force_save_frequency == 20
    assert method.thermodynamics_save_frequency == 100


@pytest.mark.unit
def test_minimize_after_run_is_not_misclassified_as_md():
    settings = selected_run(
        ['fix 1 all nve\n', 'run 10\n', 'minimize 1e-4 1e-6 100 1000\n']
    )
    assert settings == {'kind': 'minimize'}


@pytest.mark.integration
def test_log_only_run_has_serializable_workflow(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\nfix 1 all nve\ntimestep 1\nrun 10\n'
    )
    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.workflow2.method.thermodynamic_ensemble == 'NVE'
    assert archive.workflow2.method.n_steps == 10
    assert archive.workflow2.method.integration_timestep.to('fs').magnitude == approx(1)
    assert archive.m_to_dict()['workflow2']['method']['n_steps'] == 10


@pytest.mark.integration
def test_log_only_minimization_does_not_create_md_workflow(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\nfix 1 all nve\nrun 10\n'
        'minimize 1e-4 1e-6 100 1000\n'
    )
    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.workflow2 is None


@pytest.mark.integration
def test_log_only_separate_temperature_and_pressure_controls_are_serializable(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\n'
        'fix md all nve\n'
        'fix bath all temp/csld 300 300 100 12345\n'
        'fix baro all press/berendsen aniso 1 1 1000\n'
        'run 10\n'
    )
    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    method = archive.workflow2.method
    assert method.thermodynamic_ensemble == 'NPT'
    assert method.thermostat_parameters[0].thermostat_type == (
        'velocity_rescaling_langevin'
    )
    assert method.barostat_parameters[0].barostat_type == 'berendsen'
    assert method.barostat_parameters[0].coupling_type == 'anisotropic'
    assert np.diag(
        method.barostat_parameters[0].reference_pressure.to('bar').magnitude
    ) == approx([1, 1, 1])
    assert archive.m_to_dict()['workflow2']['method']['n_steps'] == 10
