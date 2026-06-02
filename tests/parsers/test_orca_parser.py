from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
)

from nomad_simulation_parsers.parsers.orca.parser import OrcaParser

LOGGER = get_logger(__name__)
DATA_DIR = Path(__file__).resolve().parents[1] / 'data' / 'orca'


def _parse_orca(filename: str) -> EntryArchive:
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse(str(DATA_DIR / filename), archive, LOGGER)
    return archive


def test_parse_file():
    archive = _parse_orca('RI_MP2_water.out')

    assert archive.data.program.name == 'ORCA'
    assert archive.data.program.version == '5.0.4'


def test_model_system_and_molecular_orbitals():
    archive = _parse_orca('RI_MP2_water.out')

    assert len(archive.data.model_system) == 1
    system = archive.data.model_system[0]
    assert system.is_representative
    assert [state.chemical_symbol for state in system.particle_states] == [
        'O',
        'H',
        'H',
    ]
    assert system.positions[0][2].to('angstrom').magnitude == pytest.approx(0.11779)
    assert system.total_charge == 0
    assert system.total_spin == 0

    assert len(archive.data.outputs) == 1
    output = archive.data.outputs[0]
    assert output.model_system_ref.m_proxy_value == '/data/model_system/0'
    assert len(output.electronic_eigenvalues) == 1

    molecular_orbitals = output.electronic_eigenvalues[0]
    assert isinstance(molecular_orbitals, MolecularOrbitals)
    assert molecular_orbitals.n_mo == 24
    assert molecular_orbitals.n_ao == 24
    assert molecular_orbitals.mo_type == 'canonical'
    assert molecular_orbitals.mo_occupations[0] == pytest.approx(2.0)
    assert molecular_orbitals.mo_occupations[-1] == pytest.approx(0.0)
    assert molecular_orbitals.mo_energies[0].to('electron_volt').magnitude == (
        pytest.approx(-559.0826)
    )


def test_molecular_orbital_coefficients():
    archive = _parse_orca('orca_orbitals.out')

    molecular_orbitals = archive.data.outputs[0].electronic_eigenvalues[0]
    assert molecular_orbitals.n_mo == 77
    assert molecular_orbitals.n_ao == 77
    assert molecular_orbitals.mo_coefficients.shape == (77, 77)
    assert molecular_orbitals.mo_coefficients[0, 0] == pytest.approx(-0.000034)
    assert molecular_orbitals.mo_coefficients[5, 9] == pytest.approx(0.995146)
    assert molecular_orbitals.mo_coefficients[-1, -1] == pytest.approx(-0.006075)
    assert molecular_orbitals.mo_occupations[28] == pytest.approx(0.0)
    assert molecular_orbitals.mo_energies[0].to('electron_volt').magnitude == (
        pytest.approx(-520.4853)
    )
