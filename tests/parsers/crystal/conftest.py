from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.crystal.parser import CrystalParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'crystal'
LOGGER = get_logger(__name__)


def parse_crystal(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    CrystalParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def single_point_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'single_point' / 'dft' / 'output.out')


@pytest.fixture(scope='module')
def geometry_optimization_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'geo_opt' / 'nio_tzvp_pbe0_opt.o')


@pytest.fixture(scope='module')
def band_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'band_structure' / 'nacl_hf' / 'NaCl.out')


@pytest.fixture(scope='module')
def dos_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'dos' / 'nacl_hf' / 'NaCl.out')


@pytest.fixture(scope='module')
def molecule_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'molecule' / 'w.out')


@pytest.fixture(scope='module')
def band_no_f25_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'band_structure' / 'no_f25_1' / 'test04_dft.out')


@pytest.fixture(scope='module')
def band_no_f25_variant_archive() -> EntryArchive:
    return parse_crystal(
        DATA_DIR / 'band_structure' / 'no_f25_2' / 'TiS2_band_structure.prop.o'
    )


@pytest.fixture(scope='module')
def geometry_optimization_nio_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'geo_opt' / 'nio_tzvp_pbe0_opt.o')


@pytest.fixture(scope='module')
def constraints_archive() -> EntryArchive:
    return parse_crystal(
        DATA_DIR / 'misc' / 'constraints' / 'ionic1_fullspin_spinfx_2.cryst.out'
    )


@pytest.fixture(scope='module')
def displacement_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'misc' / 'displacement' / 'fe50_x8_l0_re.cryst.out')


@pytest.fixture(scope='module')
def ghosts_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'misc' / 'ghosts' / 'fevo46_sngt_ti_zero.cryst.out')


@pytest.fixture(scope='module')
def nat_archive() -> EntryArchive:
    return parse_crystal(
        DATA_DIR / 'misc' / 'nat' / 'HfS2_PBE0D3_ZD_fc3_supercell-00497.o'
    )


@pytest.fixture(scope='module')
def substitution_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'misc' / 'substitution' / 'neutral.cryst.out')


@pytest.fixture(scope='module')
def nanotube_geometry_optimization_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'nanotube' / 'geo_opt' / 'test_nano05.out')


@pytest.fixture(scope='module')
def nanotube_scf_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'nanotube' / 'scf' / 'test_nano07_3.out')


@pytest.fixture(scope='module')
def single_point_forces_archive() -> EntryArchive:
    return parse_crystal(
        DATA_DIR / 'single_point' / 'forces' / 'HfS2_PBE0D3_ZD_fc3_supercell-00001.o'
    )


@pytest.fixture(scope='module')
def single_point_hf_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'single_point' / 'hf' / 'output.out')


@pytest.fixture(scope='module')
def surface_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'surface' / 'w221_sr_pbe0.cryst.out')


@pytest.fixture(scope='module')
def pbe0_archive() -> EntryArchive:
    return parse_crystal(
        DATA_DIR / 'xc_functionals' / 'pbe0' / 'ZrS2_band_structure_dos.prop.o'
    )


@pytest.fixture(scope='module')
def pbe_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'xc_functionals' / 'pbe_1' / 'supercell-00138.o')


@pytest.fixture(scope='module')
def pw91_hybrid_archive() -> EntryArchive:
    return parse_crystal(DATA_DIR / 'xc_functionals' / 'pw91_hybrid' / 'f075_l3_ph.o')


@pytest.fixture(scope='module')
def wc1lyp_archive() -> EntryArchive:
    return parse_crystal(
        DATA_DIR / 'xc_functionals' / 'wc1lyp' / 'albite_freq_intens.out'
    )
