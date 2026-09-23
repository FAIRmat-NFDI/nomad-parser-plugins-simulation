from pathlib import Path

import phonopy
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.phonopy.parser import (
    PhonopyParser,
    phonopy_obj_to_archive,
)

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'phonopy'
LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def vasp_phonopy_archive() -> EntryArchive:
    archive = EntryArchive()
    PhonopyParser().parse(str(DATA_DIR / 'vasp' / 'phonopy.yaml'), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def cp2k_phonopy_object():
    return phonopy.load(
        DATA_DIR / 'cp2k_hexagonal_noncanonical' / 'phonopy.yaml',
        force_sets_filename=DATA_DIR / 'cp2k_hexagonal_noncanonical' / 'FORCE_SETS',
    )


@pytest.fixture(scope='module')
def cp2k_phonopy_archive(cp2k_phonopy_object) -> EntryArchive:
    archive = EntryArchive()
    phonopy_obj_to_archive(cp2k_phonopy_object, archive, LOGGER)
    return archive
