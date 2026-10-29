from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_file_parser.mapping_parser import TextParser

from nomad_simulation_parsers.parsers import cp2k_parser
from nomad_simulation_parsers.parsers.cp2k.file_parser import CP2KOutParser
from nomad_simulation_parsers.parsers.cp2k.parser import CP2KOutputParser

DATA = Path(__file__).resolve().parents[2] / 'data/cp2k/single_point/si_bulk8.out'


@pytest.mark.unit
def test_source_transformations_use_raw_reader_data():
    parser = CP2KOutputParser(filepath=str(DATA), text_parser=CP2KOutParser())
    source = parser.data
    assert CP2KOutputParser.to_dict is TextParser.to_dict
    assert 'quickstep' in source
    assert 'systems' not in source
    assert 'outputs' not in source
    assert parser.get_version(source) == 'CP2K version 2.6.2'
    assert parser.get_run_type(source) == 'ENERGY_FORCE'
    system = parser.get_systems(source)[0]
    assert len(system['atoms']) == 8
    assert system['positions'].to('angstrom').magnitude[6, 2] == pytest.approx(4.073023)
    output = parser.get_outputs(source)[0]
    assert output['energy'].to('hartree').magnitude == pytest.approx(
        -31.297885372811063
    )
    assert output['forces'].to('hartree / bohr').shape == (8, 3)
    assert len(output['scf_steps']['energies_total']) == 10


@pytest.mark.unit
def test_dftb_fallback_and_missing_data():
    parser = CP2KOutputParser()
    coordinates = SimpleNamespace(
        positions=np.zeros((1, 3)),
        labels=['Si'],
        atomic_numbers={'Si': 14},
    )
    source = {
        'global': [['run_type', 'ENERGY']],
        'qs_dftb': {
            'atomic_coordinates': coordinates,
            'single_point': {},
        },
        'lattice_vectors': [[1, 0, 0]],
    }
    systems = parser.get_systems(source)
    assert systems[0]['atoms'] == [{'label': 'Si', 'number': 14}]
    assert 'lattice_vectors' not in systems[0]
    assert parser.get_outputs(source) == []
    assert parser.get_systems({'qs_dftb': {'atomic_coordinates': None}}) == []
    assert parser.get_outputs({'global': [['run_type', 'ENERGY']]}) == []


@pytest.mark.unit
def test_recognition(tmp_path):
    contents = DATA.read_text()
    parser = cp2k_parser.load()
    assert parser.is_mainfile(str(DATA), 'text/plain', contents.encode(), contents)
    other = tmp_path / 'other.out'
    other.write_text('not a CP2K output')
    assert not parser.is_mainfile(
        str(other), 'text/plain', b'not a CP2K output', 'not a CP2K output'
    )


@pytest.mark.integration
def test_single_point():
    archive = EntryArchive()
    cp2k_parser.load().parse(str(DATA), archive, get_logger(__name__))

    simulation = archive.data
    assert simulation.program.name == 'CP2K'
    assert simulation.program.version == 'CP2K version 2.6.2'
    system = simulation.model_system[0]
    assert [atom.chemical_symbol for atom in system.particle_states] == ['Si'] * 8
    assert [atom.atomic_number for atom in system.particle_states] == [14] * 8
    assert system.positions.to('angstrom').magnitude[6, 2] == pytest.approx(4.073023)
    np.testing.assert_allclose(
        system.lattice_vectors.to('angstrom').magnitude, np.eye(3) * 5.431
    )
    output = simulation.outputs[0]
    assert output.total_energies[0].value.to('hartree').magnitude == pytest.approx(
        -31.297885372811063
    )
    assert output.total_forces[0].value.shape == (8, 3)
    assert len(output.scf_steps.energies_total) == 10
    assert archive.workflow2.m_def.name == 'SinglePoint'
    assert archive.m_validate() == ([], [])


@pytest.mark.integration
def test_reuse_with_unsupported_run_type(tmp_path):
    parser = cp2k_parser.load()
    logger = get_logger(__name__)
    first = EntryArchive()
    parser.parse(str(DATA), first, logger)

    other_file = tmp_path / 'geo_opt.out'
    other_file.write_text(DATA.read_text().replace('ENERGY_FORCE', 'GEO_OPT'))
    second = EntryArchive()
    parser.parse(str(other_file), second, logger)

    assert first.workflow2.m_def.name == 'SinglePoint'
    assert second.workflow2 is None
    assert len(second.data.model_system) == 1
    assert not second.data.outputs


@pytest.mark.integration
def test_scf_iterations_with_large_orbital_table(tmp_path):
    mainfile = tmp_path / 'orbitals.out'
    orbital_rows = (
        ' MO|     18              -2.830094'
        '             -77.010781               2.000000\n' * 2000
    )
    mainfile.write_text(
        DATA.read_text().replace(
            'SCF WAVEFUNCTION OPTIMIZATION',
            f'SCF WAVEFUNCTION OPTIMIZATION\n{orbital_rows}',
            1,
        )
    )
    archive = EntryArchive()
    cp2k_parser.load().parse(str(mainfile), archive, get_logger(__name__))
    assert len(archive.data.outputs[0].scf_steps.energies_total) == 10


@pytest.mark.integration
def test_missing_quickstep_does_not_fabricate_results(tmp_path):
    mainfile = tmp_path / 'incomplete.out'
    mainfile.write_text(' CP2K| version string:    CP2K version 2024.1\n')
    archive = EntryArchive()
    cp2k_parser.load().parse(str(mainfile), archive, get_logger(__name__))
    assert archive.data.program.version == 'CP2K version 2024.1'
    assert not archive.data.model_system
    assert not archive.data.outputs
    assert archive.workflow2 is None
