from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lobster.parser import LobsterParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'lobster'
LOGGER = get_logger(__name__)


def parse_lobster(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    LobsterParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def fe_archive() -> EntryArchive:
    return parse_lobster(DATA_DIR / 'Fe' / 'lobsterout')


@pytest.fixture(scope='module')
def nacl_archive() -> EntryArchive:
    return parse_lobster(DATA_DIR / 'NaCl' / 'lobsterout')


@pytest.fixture(scope='module')
def hfv2_archive() -> EntryArchive:
    return parse_lobster(DATA_DIR / 'HfV2' / 'lobsterout')


@pytest.fixture(scope='module')
def ni_archive() -> EntryArchive:
    return parse_lobster(DATA_DIR / 'Ni' / 'lobsterout')


@pytest.fixture(scope='module')
def si_archive() -> EntryArchive:
    return parse_lobster(DATA_DIR / 'Si' / 'lobsterout.gz')
