"""Reject reduced LAMMPS units before assigning physical schema quantities."""

from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.file_parsers import get_unit
from nomad_simulation_parsers.parsers.lammps.parser import (
    LammpsArchiveWriter,
    LammpsParser,
)

LOGGER = get_logger(__name__)
DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'lammps'


@pytest.mark.integration
def test_rejects_lj_fixture_before_writing_archive():
    archive = EntryArchive()
    mainfile = DATA_DIR / '1_xyz_files' / 'log.lammps'

    with pytest.raises(ValueError, match='LAMMPS units lj are reduced'):
        LammpsParser().parse(str(mainfile), archive, LOGGER)

    assert archive.data is None
    assert archive.workflow2 is None
    assert get_unit('lj') == {}


@pytest.mark.unit
def test_writer_cannot_bypass_lj_rejection():
    archive = EntryArchive()
    mainfile = DATA_DIR / '1_xyz_files' / 'log.lammps'

    with pytest.raises(ValueError, match='LAMMPS units lj are reduced'):
        LammpsArchiveWriter().write(str(mainfile), archive, LOGGER)

    assert archive.data is None


@pytest.mark.unit
@pytest.mark.parametrize(
    ('contents', 'message'),
    [
        ('LAMMPS (test)\nrun 10\n', 'unit style is unresolved'),
        ('LAMMPS (test)\nunits real\nunits metal\n', 'Conflicting LAMMPS'),
        (
            'LAMMPS (test)\nunits real\n  Unit style    : metal\n',
            'Conflicting LAMMPS',
        ),
        (
            'LAMMPS (test)\n  Unit style    : real\n  Unit style    : metal\n',
            'Conflicting LAMMPS',
        ),
        ('LAMMPS (test)\n  Unit style    : lj\n', 'LAMMPS units lj are reduced'),
        ('LAMMPS (test)\nunits imaginary\n', 'Unsupported LAMMPS'),
        ('LAMMPS (test)\n  Unit style    : imaginary\n', 'Unsupported LAMMPS'),
    ],
)
def test_rejects_ambiguous_units_before_writing_archive(tmp_path, contents, message):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(contents)
    archive = EntryArchive()

    with pytest.raises(ValueError, match=message):
        LammpsParser().parse(str(mainfile), archive, LOGGER)

    assert archive.data is None
    assert archive.workflow2 is None


@pytest.mark.unit
@pytest.mark.parametrize(
    ('unit_style', 'expected_velocity'), [('real', 20000), ('metal', 20)]
)
def test_physical_unit_styles_keep_velocity_scales(unit_style, expected_velocity):
    velocity = 0.2 * get_unit(unit_style, 'velocity')
    assert velocity.to('meter/second').magnitude == pytest.approx(expected_velocity)


@pytest.mark.unit
@pytest.mark.parametrize('unit_style', ['real', 'metal'])
def test_same_log_runtime_units_are_used_for_conversion(tmp_path, unit_style):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(f'LAMMPS (test)\n  Unit style    : {unit_style}\n')
    writer = LammpsArchiveWriter()
    writer.mainfile = str(mainfile)

    writer._resolved_unit_style = writer._require_physical_unit_style()
    writer._configure_parsers()

    assert writer._resolved_unit_style == unit_style
    assert writer._log_parser.units == get_unit(unit_style)
    writer._log_parser.close()


@pytest.mark.integration
@pytest.mark.parametrize('unit_style', ['real', 'metal'])
def test_runtime_only_units_are_used_for_md_workflow(tmp_path, unit_style):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        f'LAMMPS (test)\n  Unit style    : {unit_style}\n'
        'fix md all nvt temp 300 300 100\n'
        'timestep 0.5\nrun 10\n'
    )
    archive = EntryArchive()

    LammpsParser().parse(str(mainfile), archive, LOGGER)

    method = archive.workflow2.method
    time_unit = 'fs' if unit_style == 'real' else 'ps'
    assert method.thermodynamic_ensemble == 'NVT'
    assert method.integration_timestep.to(time_unit).magnitude == pytest.approx(0.5)
    assert method.thermostat_parameters[0].coupling_constant.to(
        time_unit
    ).magnitude == pytest.approx(100)
    assert 'integration_timestep' in archive.m_to_dict()['workflow2']['method']


@pytest.mark.unit
def test_matching_command_and_runtime_units_are_accepted(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text('LAMMPS (test)\nunits real\n  Unit style    : real\n')
    writer = LammpsArchiveWriter()
    writer.mainfile = str(mainfile)

    assert writer._require_physical_unit_style() == 'real'


@pytest.mark.unit
def test_neighboring_input_does_not_resolve_log_units(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text('LAMMPS (test)\nrun 10\n')
    (tmp_path / 'in.lammps').write_text('units real\nrun 10\n')
    archive = EntryArchive()

    with pytest.raises(ValueError, match='unit style is unresolved'):
        LammpsParser().parse(str(mainfile), archive, LOGGER)

    assert archive.data is None
