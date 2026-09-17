from pathlib import Path

import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'h5md'
MAINFILE = DATA_DIR / 'test_traj_openmm_reduced-SOL_5frames_07-10-25.h5'


@pytest.fixture(scope='module')
def parser() -> H5MDParser:
    return H5MDParser()


@pytest.fixture(scope='module')
def h5md_archive() -> EntryArchive:
    archive = EntryArchive()
    H5MDParser().parse(str(MAINFILE), archive, get_logger(__name__))
    normalize_all(archive, logger=get_logger(__name__))
    return archive
