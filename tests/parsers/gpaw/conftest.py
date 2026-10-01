from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.gpaw.parser import GPAWParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'gpaw'
LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def parser() -> GPAWParser:
    return GPAWParser()


def parse_gpaw(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    GPAWParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def fe2_archive() -> EntryArchive:
    return parse_gpaw(DATA_DIR / 'Fe2.gpw')


@pytest.fixture(scope='module')
def h2_archive() -> EntryArchive:
    return parse_gpaw(DATA_DIR / 'H2.gpw')


@pytest.fixture(scope='module')
def hspinpol_archive() -> EntryArchive:
    return parse_gpaw(DATA_DIR / 'Hspinpol.gpw')


@pytest.fixture(scope='module')
def si_pw_archive() -> EntryArchive:
    return parse_gpaw(DATA_DIR / 'Si_pw.gpw2')


@pytest.fixture(scope='module')
def si_lcao_archive() -> EntryArchive:
    return parse_gpaw(DATA_DIR / 'Si_lcao.gpw2')
