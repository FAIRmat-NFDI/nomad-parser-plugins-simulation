from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.gaussian.parser import GaussianParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'gaussian'
LOGGER = get_logger(__name__)


def parse_gaussian(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    GaussianParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def al_scf_archive() -> EntryArchive:
    return parse_gaussian(DATA_DIR / 'Al_scf' / 'Al.out')


@pytest.fixture(scope='module')
def al_multistep_archive() -> EntryArchive:
    return parse_gaussian(DATA_DIR / 'Al_multistep' / 'm61b5.out')


@pytest.fixture(scope='module')
def cho_frequency_archive() -> EntryArchive:
    return parse_gaussian(DATA_DIR / 'CHO_freq' / 'prono.out')


@pytest.fixture(scope='module')
def no_mp_archive() -> EntryArchive:
    return parse_gaussian(DATA_DIR / 'NO_mp' / 'onno.out')
