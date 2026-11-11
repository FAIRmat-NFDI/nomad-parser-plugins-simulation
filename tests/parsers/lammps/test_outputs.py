"""Source-backed LAMMPS thermo output regressions."""

import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger
from nomad_simulations.schema_packages.outputs import (
    TrajectoryOutputs,
    WorkflowOutputs,
)

from nomad_simulation_parsers.parsers.lammps.file_parsers import LogParser
from tests.parsers.lammps.conftest import DATA_DIR, parse_lammps


@pytest.mark.unit
@pytest.mark.parametrize(
    ('relative_path', 'n_records', 'last_step'),
    [
        (
            'hexane_cyclohexane/log.hexane_cyclohexane_nvt-thermo_style_multi',
            201,
            80000,
        ),
        ('1_methyl_naphthalene/log.naph_298_396_20ns', 301, 30000),
    ],
)
def test_thermo_reader_preserves_order_and_formats(relative_path, n_records, last_step):
    reader = LogParser()
    reader.mainfile = str(DATA_DIR / relative_path)

    blocks = reader.get_thermodynamic_blocks()

    assert len(blocks) == 1
    assert len(blocks[0]) == n_records
    assert blocks[0][0]['Step'] == 0
    assert blocks[0][-1]['Step'] == last_step
    if relative_path.endswith('thermo_style_multi'):
        assert blocks[0][0]['CPU'] == 0
        assert blocks[0][1]['CPU'] == pytest.approx(1.0376)


@pytest.mark.integration
def test_hexane_thermo_outputs_keep_values_and_system_links(hexane_archive):
    outputs = hexane_archive.data.outputs
    systems = hexane_archive.data.model_system

    assert len(outputs) == len(systems) == 201
    assert systems[21].time_step == 8400
    assert outputs[0].step == 0
    assert outputs[-1].step == 80000
    assert outputs[21].model_system_ref is systems[21]
    assert outputs[21].time.to('fs').magnitude == pytest.approx(2100)
    assert outputs[1].cpu1_end.to('second').magnitude == pytest.approx(1.0376)
    assert outputs[21].total_energies[0].value.to('joule').magnitude == (
        pytest.approx(8.866891968814006e-18)
    )
    contributions = {
        contribution.name: contribution
        for contribution in outputs[21].total_energies[0].contributions
    }
    assert contributions['PotEng'].contribution_type == 'potential'
    assert contributions['E_vdwl'].contribution_type == 'van der Waals'
    assert contributions['E_coul'].contribution_type == 'coulomb'
    assert contributions['E_bond'].value.to('joule').magnitude == pytest.approx(
        1.917461089369534e-18
    )


@pytest.mark.integration
def test_methyl_thermo_outputs_do_not_invent_unsampled_systems(
    methyl_naphthalene_archive,
):
    outputs = methyl_naphthalene_archive.data.outputs
    systems = methyl_naphthalene_archive.data.model_system

    assert len(outputs) == 301
    assert len(systems) == 4
    assert [system.time_step for system in systems] == [0, 10000, 20000, 30000]
    assert outputs[0].step == 0
    assert outputs[-1].step == 30000
    assert [
        output.step for output in outputs if output.model_system_ref is not None
    ] == [0, 10000, 20000, 30000]
    assert all(
        output.model_system_ref in systems
        for output in outputs
        if output.model_system_ref is not None
    )
    assert outputs[0].pressures[0].value.to('bar').magnitude == pytest.approx(
        17120.9764425
    )
    assert outputs[0].volumes[0].value.to('angstrom**3').magnitude == (
        pytest.approx(19152.149)
    )
    assert outputs[0].mass_densities[0].value.to('g/cm**3').magnitude == (
        pytest.approx(0.66577573)
    )

    restored = EntryArchive.m_from_dict(methyl_naphthalene_archive.m_to_dict())
    restored.metadata = EntryMetadata()
    normalize_all(restored, logger=get_logger(__name__))
    assert [
        output.step
        for output in restored.data.outputs
        if output.model_system_ref is not None
    ] == [0, 10000, 20000, 30000]


@pytest.mark.integration
def test_dcd_frame_ordinals_are_not_assumed_to_be_lammps_steps(
    methane_dcd_archive,
):
    assert len(methane_dcd_archive.data.outputs) == 201
    assert len(methane_dcd_archive.data.model_system) == 201
    assert all(
        system.time_step is None for system in methane_dcd_archive.data.model_system
    )
    assert all(
        output.model_system_ref is None for output in methane_dcd_archive.data.outputs
    )


@pytest.mark.integration
def test_log_only_thermo_outputs_round_trip_and_normalize(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\ntimestep 1\nfix int all nve\nrun 10\n'
        'Step Temp PotEng KinEng TotEng Press\n'
        '0 300 -5 2 -3 1\n'
        '5 310 -4 3 -1 2\n'
        '10 320 -3 4 1 3\n'
        'Loop time of 1 on 1 procs for 10 steps\n'
    )
    archive = parse_lammps(mainfile)

    assert len(archive.data.model_system) == 0
    assert [output.step for output in archive.data.outputs] == [0, 5, 10]
    assert [output.time.to('fs').magnitude for output in archive.data.outputs] == (
        pytest.approx([0, 5, 10])
    )
    output = archive.data.outputs[1]
    assert output.model_system_ref is None
    assert output.temperatures[0].value.to('kelvin').magnitude == 310
    assert output.pressures[0].value.to('pascal').magnitude == pytest.approx(202650)
    assert output.total_energies[0].value.to('joule').magnitude == pytest.approx(
        -6.947695345147254e-21
    )
    assert {
        contribution.contribution_type
        for contribution in output.total_energies[0].contributions
    } == {'potential', 'kinetic'}

    restored = EntryArchive.m_from_dict(archive.m_to_dict())
    restored.metadata = EntryMetadata()
    normalize_all(restored, logger=get_logger(__name__))
    assert [output.step for output in restored.data.outputs] == [0, 5, 10]


@pytest.mark.integration
def test_multiple_thermo_blocks_keep_only_selected_final_run(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\nfix int all nve\n'
        'run 10\nStep Temp\n0 300\n10 301\n'
        'Loop time of 1 on 1 procs for 10 steps\n'
        'run 20\nStep Temp\n10 302\n20 303\n'
        'Loop time of 1 on 1 procs for 20 steps\n'
    )

    archive = parse_lammps(mainfile)

    assert [output.step for output in archive.data.outputs] == [10, 20]
    assert archive.data.outputs[0].temperatures[0].value.magnitude == 302
    assert all(output.time is None for output in archive.data.outputs)


@pytest.mark.integration
def test_changed_thermo_style_within_run_preserves_all_rows(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\nfix int all nve\nrun 10\n'
        'Step Temp\n0 300\n5 301\n'
        'Step Temp PotEng\n10 302 -5\n'
        'Loop time of 1 on 1 procs for 10 steps\n'
    )

    archive = parse_lammps(mainfile)

    assert [output.step for output in archive.data.outputs] == [0, 5, 10]
    assert not archive.data.outputs[0].potential_energies
    assert archive.data.outputs[2].potential_energies[0].value.to(
        'electron_volt'
    ).magnitude == pytest.approx(-5)


@pytest.mark.integration
def test_reset_steps_do_not_infer_continuous_physical_time(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\ntimestep 1\nfix int all nve\nrun 10\n'
        'Step Temp\n0 300\n10 301\n0 302\n'
        'Loop time of 1 on 1 procs for 10 steps\n'
    )

    archive = parse_lammps(mainfile)

    assert [output.step for output in archive.data.outputs] == [0, 10, 0]
    assert all(output.time is None for output in archive.data.outputs)
    assert all(output.model_system_ref is None for output in archive.data.outputs)


@pytest.mark.integration
def test_explicit_thermo_time_survives_reset_steps(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\nfix int all nve\nrun 10\n'
        'Step Time Temp\n10 0.5 300\n0 0.6 301\n'
        'Loop time of 1 on 1 procs for 10 steps\n'
    )

    archive = parse_lammps(mainfile)

    assert [output.step for output in archive.data.outputs] == [10, 0]
    times = [output.time.to('picosecond').magnitude for output in archive.data.outputs]
    assert times == pytest.approx([0.5, 0.6])


@pytest.mark.integration
def test_minimization_thermo_rows_are_not_labeled_as_md_trajectory(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\nminimize 1e-4 1e-6 10 100\n'
        'Step PotEng E_bond TotEng Enthalpy\n'
        '0 -5 -2 -4 -3\n2 -6 -3 -5 -4\n'
        'Loop time of 1 on 1 procs for 2 steps\n'
    )

    archive = parse_lammps(mainfile)

    assert archive.workflow2 is None
    assert [output.step for output in archive.data.outputs] == [0, 2]
    assert all(isinstance(output, WorkflowOutputs) for output in archive.data.outputs)
    assert not any(
        isinstance(output, TrajectoryOutputs) for output in archive.data.outputs
    )
    assert archive.data.outputs[1].total_energies[0].value.to(
        'electron_volt'
    ).magnitude == pytest.approx(-5)
    assert archive.data.outputs[1].enthalpies[0].value.to(
        'electron_volt'
    ).magnitude == pytest.approx(-4)


@pytest.mark.integration
def test_energy_components_are_retained_without_total_energy(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\nfix int all nve\nrun 1\n'
        'Step PotEng E_bond E_angle\n0 -5 1 2\n'
        'Loop time of 1 on 1 procs for 1 steps\n'
    )

    archive = parse_lammps(mainfile)
    output = archive.data.outputs[0]

    assert not output.total_energies
    assert [energy.name for energy in output.potential_energies] == [
        'PotEng',
        'E_bond',
        'E_angle',
    ]
    assert output.potential_energies[1].value.to('electron_volt').magnitude == 1


@pytest.mark.integration
def test_source_step_matched_dump_forces_are_stored_with_physical_units(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\ntimestep 0.1\nfix int all nve\n'
        'dump 1 all custom 10 trajectory.dump id type x y z fx fy fz\n'
        'run 10\nStep Temp\n0 300\n10 301\n'
        'Loop time of 1 on 1 procs for 10 steps\n'
    )
    (tmp_path / 'trajectory.dump').write_text(
        'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n'
        'ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n'
        'ITEM: ATOMS id type x y z fx fy fz\n1 1 1 2 3 1 2 3\n'
        'ITEM: TIMESTEP\n10\nITEM: NUMBER OF ATOMS\n1\n'
        'ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n'
        'ITEM: ATOMS id type x y z fx fy fz\n1 1 2 2 3 4 5 6\n'
    )

    archive = parse_lammps(mainfile)

    assert len(archive.data.model_system) == 2
    assert [output.step for output in archive.data.outputs] == [0, 10]
    assert all(output.model_system_ref is not None for output in archive.data.outputs)
    assert archive.data.outputs[0].total_forces[0].value.shape == (1, 3)
    assert archive.data.outputs[1].total_forces[0].value[0][0].to(
        'newton'
    ).magnitude == pytest.approx(4 * 1.602176634e-9)
