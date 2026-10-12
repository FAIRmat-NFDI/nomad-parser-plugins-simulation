from types import SimpleNamespace

import numpy as np
import phonopy
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.phonopy import parser as parser_module
from nomad_simulation_parsers.parsers.phonopy.calculator import PhononProperties
from nomad_simulation_parsers.parsers.phonopy.parser import (
    create_system,
    phonopy_obj_to_archive,
)
from tests.parsers.phonopy.conftest import LARGE_DATA_DIR


@pytest.mark.unit
def test_create_system_populates_structure_and_supercell():
    system = create_system(
        np.eye(3), ['Si', 'Si'], np.zeros((2, 3)), np.diag([2, 2, 2])
    )

    assert [state.chemical_symbol for state in system.particle_states] == ['Si', 'Si']
    np.testing.assert_allclose(system.positions.to('angstrom').magnitude, 0)
    np.testing.assert_allclose(
        system.representations[0].lattice_vectors.to('angstrom').magnitude,
        np.eye(3),
    )
    np.testing.assert_array_equal(
        system.representations[0].supercell_matrix, np.diag([2, 2, 2])
    )


class FakeAtoms:
    cell = np.eye(3)
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    symbols = ['Si', 'Si']


class FakePhonopyObject:
    unitcell = FakeAtoms()
    supercell = FakeAtoms()
    supercell_matrix = np.diag([1, 1, 1])
    displacements = np.array([[[0.0, 0.0, 0.0], [0.01, 0.0, 0.0]]])
    force_constants = np.zeros((2, 2, 3, 3))
    symmetry = SimpleNamespace(tolerance=1e-5)
    calculator = 'vasp'
    nac_params = None


@pytest.mark.unit
def test_phonopy_obj_to_archive_creates_simulation_systems(monkeypatch):
    archive = EntryArchive()
    properties = SimpleNamespace(
        mesh=[2, 2, 2],
        frequencies=np.array([-1.0, 1.0]),
        phonopy_obj=FakePhonopyObject(),
    )
    monkeypatch.setattr(
        parser_module, 'PhononProperties', lambda *args, **kwargs: properties
    )
    monkeypatch.setattr(parser_module, 'get_bandstructures', lambda _: [])
    monkeypatch.setattr(
        parser_module,
        'get_dos',
        lambda _: [dict(frequencies=np.array([0.0, 1.0]), dos=np.array([1.0, 1.0]))],
    )
    monkeypatch.setattr(
        parser_module,
        'get_thermodynamic_properties',
        lambda _: [dict(temperature=300.0, free_energy=2.0, heat_capacity=3.0)],
    )

    result = phonopy_obj_to_archive(FakePhonopyObject(), archive, get_logger(__name__))

    assert result is archive
    assert archive.data.program.name == 'Phonopy'
    assert archive.data.program.version
    assert len(archive.data.model_system) == 2
    assert archive.data.model_system[0].positions.shape == (2, 3)
    assert archive.data.model_system[1].positions.shape == (2, 3)
    assert len(archive.data.model_method) == 1
    assert archive.data.model_method[0].name == 'harmonic lattice dynamics'
    assert len(archive.data.outputs) == 1
    assert archive.data.outputs[0].model_method_ref is archive.data.model_method[0]
    assert archive.data.outputs[0].force_constants[0].value.shape == (2, 2, 3, 3)
    assert len(archive.data.outputs[0].phonon_dos) == 1
    assert archive.data.outputs[0].n_imaginary_frequencies == 1
    assert len(archive.data.outputs[0].vibrational_free_energies) == 1
    assert len(archive.data.outputs[0].vibrational_heat_capacities) == 1
    assert archive.m_validate() == ([], [])


@pytest.mark.large_fixture
def test_phonopy_fixture_loads_force_constants():
    phonopy_obj = phonopy.load(
        LARGE_DATA_DIR / 'vasp' / 'phonopy.yaml',
        force_constants_filename=LARGE_DATA_DIR / 'vasp' / 'force_constants.hdf5',
    )

    assert phonopy_obj.force_constants.shape == (48, 384, 3, 3)


@pytest.mark.large_fixture
def test_vasp_band_yaml_loads_unit_cell():
    phonopy_obj = phonopy.load(LARGE_DATA_DIR / 'vasp' / 'band.yaml')

    assert len(phonopy_obj.unitcell) == 48
    np.testing.assert_array_equal(phonopy_obj.supercell_matrix, np.diag([4, 4, 4]))


@pytest.mark.integration
def test_cp2k_fixture_band_segments_for_archive_writer(
    cp2k_hexagonal_noncanonical_phonopy_object,
):
    properties = PhononProperties(
        cp2k_hexagonal_noncanonical_phonopy_object,
        get_logger(__name__),
        k_mesh=2,
    )
    frequencies, bands, labels = properties.get_bandstructure()

    assert labels.tolist() == [
        ['Γ', 'M'],
        ['M', 'K'],
        ['K', 'Γ'],
        ['Γ', 'A'],
        ['A', 'L'],
        ['L', 'H'],
        ['H', 'A'],
        ['L', 'M'],
        ['K', 'H'],
    ]
    assert frequencies.shape == (9, 100, 54)
    assert bands.shape == (9, 100, 3)
