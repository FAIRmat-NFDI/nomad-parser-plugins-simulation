import numpy as np
import phonopy
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from tests.parsers.phonopy.conftest import DATA_DIR
from nomad_simulation_parsers.parsers.phonopy.calculator import PhononProperties
from nomad_simulation_parsers.parsers.phonopy.parser import (
    create_system,
    phonopy_obj_to_archive,
)


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


@pytest.mark.unit
def test_phonopy_obj_to_archive_creates_simulation_systems():
    archive = EntryArchive()

    result = phonopy_obj_to_archive(
        FakePhonopyObject(), archive, get_logger(__name__)
    )

    assert result is archive
    assert archive.data.program.name == 'Phonopy'
    assert archive.data.program.version
    assert len(archive.data.model_system) == 2
    assert archive.data.model_system[0].positions.shape == (2, 3)
    assert archive.data.model_system[1].positions.shape == (2, 3)
    assert archive.m_validate() == ([], [])


@pytest.mark.unit
def test_phonopy_fixture_loads_force_constants():
    phonopy_obj = phonopy.load(
        DATA_DIR / 'vasp' / 'phonopy.yaml',
        force_constants_filename=DATA_DIR / 'vasp' / 'force_constants.hdf5',
    )

    assert phonopy_obj.force_constants.shape == (48, 384, 3, 3)


@pytest.mark.integration
def test_cp2k_fixture_band_segments_for_archive_writer(cp2k_phonopy_object):
    properties = PhononProperties(cp2k_phonopy_object, get_logger(__name__), k_mesh=2)
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
