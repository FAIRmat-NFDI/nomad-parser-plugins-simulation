import numpy as np
import phonopy
import pytest
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.phonopy.calculator import PhononProperties
from tests.parsers.phonopy.conftest import LARGE_DATA_DIR


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
