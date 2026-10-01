import pytest

from nomad_simulation_parsers.parsers import gromacs_parser


@pytest.mark.unit
def test_recognizes_gromacs_mdrun_log(tmp_path):
    contents = 'gmx mdrun, VERSION 2025.0\nInput Parameters:\n'
    mainfile = tmp_path / 'md.log'
    mainfile.write_text(contents)

    parser = gromacs_parser.load()

    assert parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)


@pytest.mark.unit
@pytest.mark.parametrize(
    'contents',
    [
        'GROMACS trajectory data\n',
        'gmx mdrun, VERSION 2025.0\nno input parameters\n',
        'This is not a simulation log\n',
    ],
)
def test_rejects_non_gromacs_mainfiles(tmp_path, contents):
    mainfile = tmp_path / 'output'
    mainfile.write_text(contents)
    parser = gromacs_parser.load()

    assert not parser.is_mainfile(
        str(mainfile), 'text/plain', contents.encode(), contents
    )
