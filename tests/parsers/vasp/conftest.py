from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ServerContext
from nomad.files import StagingUploadFiles
from nomad.processing import Upload
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'vasp'
LOGGER = get_logger(__name__)


def parse_vasp(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    VASPParser().parse(str(mainfile), archive, LOGGER)
    return archive


SILICON_GW_ARCHIVE = parse_vasp(DATA_DIR / 'Si_GW' / 'vasprun.xml')


@pytest.fixture(scope='module')
def agac_relax_outcar_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'AgAc_relax' / 'OUTCAR')


@pytest.fixture(scope='module')
def agac_relax_vasprun_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'AgAc_relax' / 'vasprun.xml.relax')


@pytest.fixture(scope='module')
def mg_static_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'Mg_bands' / 'vasprun.xml.static')


@pytest.fixture(scope='module')
def mg_bands_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'Mg_bands' / 'vasprun.xml.bands')


@pytest.fixture(scope='module')
def silicon_dos_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'Si_dos_band' / 'dos_si.xml')


@pytest.fixture(scope='module')
def silicon_band_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'Si_dos_band' / 'band_si.xml')


@pytest.fixture(scope='module')
def silicon_gw_archive() -> EntryArchive:
    return SILICON_GW_ARCHIVE


@pytest.fixture(scope='module')
def silicon_gw_kspace_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'Si_GW' / 'vasprun.xml')


@pytest.fixture(scope='module')
def gamma_outcar_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'gamma' / 'OUTCAR')


@pytest.fixture
def gamma_outcar_kmesh_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'gamma' / 'OUTCAR')


@pytest.fixture(scope='module')
def al_n_vasprun_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'alternative_pseudopotentials' / 'AlN' / 'vasprun.xml')


@pytest.fixture(scope='module')
def al_n_outcar_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'alternative_pseudopotentials' / 'AlN' / 'OUTCAR')


@pytest.fixture(scope='module')
def f_gw_outcar_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'alternative_pseudopotentials' / 'F_GW' / 'OUTCAR')


@pytest.fixture(scope='module')
def hybrid_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'hybrid_vasprun.xml.gz')


@pytest.fixture(scope='module')
def metagga_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'hle17_vasprun.xml.gz')


@pytest.fixture(scope='module')
def malformed_time_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'malformed_time' / 'vasprun.xml')


@pytest.fixture(scope='module')
def broken_vasprun_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'vasprun.xml.broken')


@pytest.fixture(scope='module')
def dftu_multi_parameter_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'dftu' / 'multi_parameter' / 'vasprun.xml')


@pytest.fixture(scope='module')
def dftu_single_parameter_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'dftu' / 'single_parameter' / 'vasprun.xml')


@pytest.fixture(scope='module')
def dftu_multi_parameter_no_incar_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'dftu' / 'multi_parameter_no_incar' / 'OUTCAR')


@pytest.fixture(scope='module')
def boolean_dotted_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'booleans' / 'dotted' / 'OUTCAR')


@pytest.fixture(scope='module')
def boolean_single_char_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'booleans' / 'single_char' / 'OUTCAR')


@pytest.fixture(scope='module')
def boolean_mixed_case_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'booleans' / 'mixed_case' / 'OUTCAR')


@pytest.fixture(scope='module')
def boolean_lowercase_archive() -> EntryArchive:
    return parse_vasp(DATA_DIR / 'booleans' / 'lowercase' / 'OUTCAR')


@pytest.fixture
def chgcar_archive() -> EntryArchive:
    upload_id = 'vasp_chgcar_test_upload'
    StagingUploadFiles(upload_id=upload_id, create=True)
    upload = Upload(upload_id=upload_id)
    context = ServerContext(upload=upload)
    archive = EntryArchive(
        m_context=context,
        metadata=EntryMetadata(upload_id=upload_id, entry_id='vasp_chgcar_test'),
    )
    VASPParser().parse(str(DATA_DIR / 'with_chgcar' / 'OUTCAR'), archive, LOGGER)
    return archive
