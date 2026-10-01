from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.ams.parser import AMSParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'ams'
LOGGER = get_logger(__name__)


def parse_ams(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    AMSParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def scf_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'scf' / 'phenylrSmall-metagga.out')


@pytest.fixture(scope='module')
def band_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'band_pbe_GF' / 'band_pbe_GF.out')


@pytest.fixture(scope='module')
def dos_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'dos' / 'NiO-dos.out')


@pytest.fixture(scope='module')
def restricted_dos_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'dos' / 'NiO-dos-restricted.out')


@pytest.fixture(scope='module')
def geometry_optimization_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'go' / 'phenylrSmall-geoopt.out')


@pytest.fixture(scope='module')
def large_geometry_optimization_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'go' / 'EDUSIF.out')


@pytest.fixture(scope='module')
def adf_archive() -> EntryArchive:
    return parse_ams(DATA_DIR / 'adf_SP' / 'adf_SP.out')
