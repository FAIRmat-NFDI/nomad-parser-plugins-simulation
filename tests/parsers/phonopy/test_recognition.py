import pytest

from nomad_simulation_parsers.parsers import phonopy_parser


@pytest.mark.unit
def test_recognizes_phonopy_yaml(tmp_path):
    mainfile = tmp_path / 'phonopy.yaml'
    mainfile.write_text('phonopy:\n  version: 3.0\n')
    parser = phonopy_parser.load()

    assert parser.is_mainfile(str(mainfile), 'text/plain', b'', '')


@pytest.mark.unit
@pytest.mark.parametrize('name', ['phonopy.yml', 'POSCAR', 'phonopy.yaml.bak'])
def test_rejects_non_phonopy_mainfile_names(tmp_path, name):
    mainfile = tmp_path / name
    mainfile.write_text('phonopy:\n  version: 3.0\n')
    parser = phonopy_parser.load()

    assert not parser.is_mainfile(str(mainfile), 'text/plain', b'', '')
