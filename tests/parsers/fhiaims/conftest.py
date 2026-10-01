from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'fhiaims'
LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def parser() -> FHIAimsParser:
    return FHIAimsParser()


def parse_fhiaims(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    parser = FHIAimsParser()
    parser.parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def si_geomopt_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'Si_geomopt' / 'out.out')


@pytest.fixture(scope='module')
def cl_na_dos_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'ClNa_dos' / 'ClNa_dos.out')


@pytest.fixture(scope='module')
def fe_scf_spinpol_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'Fe_scf_spinpol' / 'out.out')


@pytest.fixture(scope='module')
def fe_band_spinpol_archive() -> EntryArchive:
    return parse_fhiaims(
        DATA_DIR / 'Fe_band_spinpol' / 'Fe_band_structure_dos_spin.out'
    )


@pytest.fixture(scope='module')
def ho_md_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'HO_md' / 'H2O_periodic_MD.out')


@pytest.fixture(scope='module')
def gaas_hse06_soc_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'GaAs_HSE06+SOC' / 'aims.out')


@pytest.fixture(scope='module')
def ceo2_dftu_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'CeO2_dftu' / 'aims_CC.out')


@pytest.fixture(scope='module')
def he_gw_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'He_gw' / 'He_scGW_ontop_PBE.out')


@pytest.fixture(scope='module')
def silicon_v071914_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'Si_band_dos_v071914_7' / 'aims_CC.out')


@pytest.fixture(scope='module')
def silicon_v171221_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'Si_band_dos_v171221_1' / 'aims_CC.out')


@pytest.fixture(scope='module')
def chn_gw_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'CHN_gw' / 'output.out')


@pytest.fixture(scope='module')
def si_gw_bands_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'Si_pbe_vs_gw_bands' / 'aims.out')


@pytest.fixture(scope='module')
def native_tight_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'native_tiers' / 'tight' / 'aims.out')


@pytest.fixture(scope='module')
def native_intermediate_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'native_tiers' / 'intermediate' / 'aims.out')


@pytest.fixture(scope='module')
def native_light_spd_archive() -> EntryArchive:
    return parse_fhiaims(DATA_DIR / 'native_tiers' / 'light_spd' / 'aims.out')
