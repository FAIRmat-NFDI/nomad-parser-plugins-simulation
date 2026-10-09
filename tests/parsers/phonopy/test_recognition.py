import numpy as np
import pytest

from nomad_simulation_parsers.parsers import phonopy_parser
from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyYamlParser


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


@pytest.mark.unit
def test_phonopy_yaml_parser_load_file_and_to_dict(tmp_path):
    mainfile = tmp_path / 'phonopy.yaml'
    mainfile.write_text('phonopy:\n  version: 3.0\nunit_cell:\n  lattice: []\n')

    parser = PhonopyYamlParser(filepath=str(mainfile))

    assert parser.load_file() is None
    assert parser.to_dict()['unit_cell']['lattice'] == []
    assert parser.data == parser.to_dict()


@pytest.mark.unit
def test_imaginary_frequencies_initializes_dos_mesh():
    class Properties:
        def get_dos(self):
            self.frequencies = np.array([[-1.0, 0.5, 2.0]])
            return np.array([]), np.array([])

    parser = PhonopyYamlParser(data_object=object())
    parser._phonon_properties = Properties()

    assert parser.get_imaginary_frequencies() == 1


@pytest.mark.unit
def test_vibrational_thermodynamics_returns_cached_key():
    parser = PhonopyYamlParser(data_object=object())
    cached = {
        'vibrational_free_energies': [{'value': np.array([1.0])}],
        'vibrational_heat_capacities': [{'value': np.array([2.0])}],
    }
    parser._thermodynamic_properties = cached

    assert (
        parser.get_vibrational_thermodynamics(key='vibrational_heat_capacities')
        is cached['vibrational_heat_capacities']
    )
