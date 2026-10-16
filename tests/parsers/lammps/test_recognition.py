import pytest

from nomad_simulation_parsers.parsers import lammps_parser


@pytest.mark.unit
def test_recognizes_lammps_mainfile(tmp_path):
    contents = 'LAMMPS (12 Dec 2018)\n'
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(contents)

    parser = lammps_parser.load()

    assert parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)


@pytest.mark.unit
@pytest.mark.parametrize(
    'contents',
    [
        'LAMMPS input script without a version header\n',
        'This is not a LAMMPS output\n',
        'GROMACS version: 2025.0\n',
    ],
)
def test_rejects_non_lammps_mainfiles(tmp_path, contents):
    mainfile = tmp_path / 'output'
    mainfile.write_text(contents)
    parser = lammps_parser.load()

    assert not parser.is_mainfile(
        str(mainfile), 'text/plain', contents.encode(), contents
    )
