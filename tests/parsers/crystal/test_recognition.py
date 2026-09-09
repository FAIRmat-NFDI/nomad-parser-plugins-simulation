import pytest

from nomad_simulation_parsers.parsers import crystal_parser


@pytest.mark.unit
def test_recognizes_crystal_header(tmp_path):
    contents = '\n * CRYSTAL14 *\n * version : 1\n'
    mainfile = tmp_path / 'calculation.out'
    mainfile.write_text(contents)
    parser = crystal_parser.load()

    assert (
        parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
        is True
    )


@pytest.mark.unit
def test_rejects_non_crystal_header(tmp_path):
    contents = 'not a CRYSTAL output\n'
    mainfile = tmp_path / 'calculation.out'
    mainfile.write_text(contents)
    parser = crystal_parser.load()

    assert not parser.is_mainfile(
        str(mainfile), 'text/plain', contents.encode(), contents
    )
