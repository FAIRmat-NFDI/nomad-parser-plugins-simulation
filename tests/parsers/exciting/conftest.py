from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'exciting'
LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def parser() -> ExcitingParser:
    return ExcitingParser()


def parse_exciting(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    ExcitingParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def c_minimal_archive() -> EntryArchive:
    return parse_exciting(DATA_DIR / 'C_minimal' / 'INFO.OUT')


@pytest.fixture(scope='module')
def c_gs_archive() -> EntryArchive:
    return parse_exciting(DATA_DIR / 'C_gs' / 'INFO.OUT')


@pytest.fixture(scope='module')
def ga_o_sodium_archive() -> EntryArchive:
    return parse_exciting(DATA_DIR / 'GaO_sodium' / 'INFO.OUT')


@pytest.fixture(scope='module')
def ga_o_strucopt_archive() -> EntryArchive:
    return parse_exciting(DATA_DIR / 'GaO_strucopt' / 'INFO.OUT')


@pytest.fixture(scope='module')
def ce_o_dos_archive() -> EntryArchive:
    return parse_exciting(DATA_DIR / 'CeO_dos' / 'INFO.OUT')


@pytest.fixture(scope='module')
def pb_i_hybrids_archive() -> EntryArchive:
    return parse_exciting(DATA_DIR / 'PbI_hybrids' / 'INFO.OUT')


EXCITING_EXAMPLE_FIXTURES = [
    'c_minimal_archive',
    'ga_o_sodium_archive',
    'ga_o_strucopt_archive',
    'pb_i_hybrids_archive',
]


@pytest.fixture(
    scope='module',
    params=EXCITING_EXAMPLE_FIXTURES,
    ids=EXCITING_EXAMPLE_FIXTURES,
)
def exciting_example_archive(request) -> EntryArchive:
    return request.getfixturevalue(request.param)
