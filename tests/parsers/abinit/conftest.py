from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'abinit'
LOGGER = get_logger(__name__)


def parse_abinit(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    AbinitParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def fe_archive() -> EntryArchive:
    return parse_abinit(DATA_DIR / 'Fe' / 'Fe.out')


@pytest.fixture(scope='module')
def h2_archive() -> EntryArchive:
    return parse_abinit(DATA_DIR / 'H2' / 'H2.out')


@pytest.fixture(scope='module')
def si_archive() -> EntryArchive:
    return parse_abinit(DATA_DIR / 'Si' / 'Si.out')


@pytest.fixture(scope='module')
def zro2_gw_archive() -> EntryArchive:
    return parse_abinit(DATA_DIR / 'ZrO2_GW' / 'A1.abo')
