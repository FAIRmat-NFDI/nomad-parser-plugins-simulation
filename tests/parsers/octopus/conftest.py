from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.octopus.parser import OctopusParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'octopus'
LOGGER = get_logger(__name__)


def parse_octopus(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    OctopusParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def si_scf_archive() -> EntryArchive:
    return parse_octopus(DATA_DIR / 'Si_scf' / 'stdout.txt')


@pytest.fixture(scope='module')
def fe_spinpol_archive() -> EntryArchive:
    return parse_octopus(DATA_DIR / 'Fe_spinpol' / 'stdout.txt')
