from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.yambo.parser import YamboParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'yambo'
LOGGER = get_logger(__name__)


def parse_yambo(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    YamboParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def hbn_archive() -> EntryArchive:
    return parse_yambo(DATA_DIR / 'hBN' / 'r-10b_1Ry_HF_and_locXC_gw0_em1d_ppa')


@pytest.fixture(scope='module')
def lif_archive() -> EntryArchive:
    return parse_yambo(DATA_DIR / 'LiF' / 'r-02_QP_PPA_em1d_ppa_HF_and_locXC_gw0')


@pytest.fixture(scope='module')
def gasb_archive() -> EntryArchive:
    return parse_yambo(DATA_DIR / 'GaSb' / 'r-02_GW_em1d_ppa_HF_and_locXC_gw0')


@pytest.fixture(scope='module')
def aluminum_archive() -> EntryArchive:
    return parse_yambo(DATA_DIR / 'Aluminum' / 'r-01_Lifetimes_em1d_life')


@pytest.fixture(scope='module')
def ch4_archive() -> EntryArchive:
    return parse_yambo(DATA_DIR / 'CH4_db_minimal' / 'r_setup')
