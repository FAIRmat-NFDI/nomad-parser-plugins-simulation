from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.wannier90.parser import Wannier90Parser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'wannier90' / 'lco_mlwf'
LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def lco_mlwf_archive() -> EntryArchive:
    archive = EntryArchive()
    Wannier90Parser().parse(str(DATA_DIR / 'lco.wout'), archive, LOGGER)
    return archive
