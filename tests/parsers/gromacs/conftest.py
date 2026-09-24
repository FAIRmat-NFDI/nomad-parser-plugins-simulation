from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.gromacs.parser import GromacsParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'gromacs'
LOGGER = get_logger(__name__)


def parse_gromacs(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    GromacsParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='class')
def water_archive() -> EntryArchive:
    return parse_gromacs(DATA_DIR / 'water' / 'reference_s.log')


@pytest.fixture(scope='class')
def protein_small_archive() -> EntryArchive:
    return parse_gromacs(DATA_DIR / 'protein_small' / 'md.log')


@pytest.fixture(scope='class')
def fep_archive() -> EntryArchive:
    return parse_gromacs(
        DATA_DIR
        / 'free_energy_calculations'
        / 'alchemical_transformation_single_run'
        / 'fep_run-7.log'
    )


@pytest.fixture(scope='class')
def fe_test_archive() -> EntryArchive:
    return parse_gromacs(DATA_DIR / 'fe_test' / 'md.log')


@pytest.fixture(scope='class')
def polymer_melt_archive() -> EntryArchive:
    return parse_gromacs(DATA_DIR / 'polymer_melt' / 'step4.0_minimization.log')


@pytest.fixture(scope='class')
def integrator_bd_archive() -> EntryArchive:
    return parse_gromacs(DATA_DIR / 'water_AA_ENUM_tests' / 'integrator-bd' / 'md.log')


@pytest.fixture(scope='class')
def integrator_mdvv_archive() -> EntryArchive:
    return parse_gromacs(
        DATA_DIR / 'water_AA_ENUM_tests' / 'integrator-mdvv' / 'md.log'
    )


@pytest.fixture(scope='class')
def integrator_sd_archive() -> EntryArchive:
    return parse_gromacs(DATA_DIR / 'water_AA_ENUM_tests' / 'integrator-sd' / 'md.log')


@pytest.fixture(scope='class')
def integrator_vrescale_archive() -> EntryArchive:
    return parse_gromacs(
        DATA_DIR
        / 'water_AA_ENUM_tests'
        / 'integrator-md'
        / 'thermostat-vrescale'
        / 'md.log'
    )


@pytest.fixture(scope='class')
def integrator_nosehoover_archive() -> EntryArchive:
    return parse_gromacs(
        DATA_DIR
        / 'water_AA_ENUM_tests'
        / 'integrator-md'
        / 'thermostat-nosehoover_barostat-parrinellorahman'
        / 'md.log'
    )
