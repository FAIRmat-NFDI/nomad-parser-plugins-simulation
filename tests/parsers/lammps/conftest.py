import os
import tempfile
from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'lammps'
LOGGER = get_logger(__name__)


def parse_lammps(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='class')
def methyl_naphthalene_archive() -> EntryArchive:
    return parse_lammps(DATA_DIR / '1_methyl_naphthalene' / 'log.1_methyl_naphthalene')


@pytest.fixture(scope='class')
def xyz_archive() -> EntryArchive:
    return parse_lammps(DATA_DIR / '1_xyz_files' / 'log.lammps')


@pytest.fixture
def tmp_dir():
    parent_directory = '.volumes'
    os.makedirs(parent_directory, exist_ok=True)
    directory = tempfile.TemporaryDirectory(dir=parent_directory, prefix='test_tmp')
    yield directory.name
    directory.cleanup()


@pytest.fixture(scope='class')
def hexane_archive() -> EntryArchive:
    return parse_lammps(DATA_DIR / 'hexane_cyclohexane' / 'log.hexane_cyclohexane_nvt')


@pytest.fixture(scope='class')
def methane_dcd_archive() -> EntryArchive:
    return parse_lammps(
        DATA_DIR / 'methane_dcd' / 'log.methane_nvt_traj_dcd_thermo_style_custom'
    )


@pytest.fixture(scope='class')
def methane_xyz_archive() -> EntryArchive:
    return parse_lammps(
        DATA_DIR / 'methane_xyz' / 'log.methane_nvt_traj_xyz_thermo_style_custom'
    )


@pytest.fixture(scope='class')
def polymer_melt_minimization_archive() -> EntryArchive:
    return parse_lammps(
        DATA_DIR / 'polymer_melt' / 'Emin' / 'log.step4.0_minimization'
    )
