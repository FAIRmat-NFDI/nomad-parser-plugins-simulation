from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoParser,
)

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'quantumespresso'
LOGGER = get_logger(__name__)


def parse_quantumespresso(mainfile: Path) -> EntryArchive:
    archive = EntryArchive()
    QuantumEspressoParser().parse(str(mainfile), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def pwscf_archive() -> EntryArchive:
    return parse_quantumespresso(DATA_DIR / 'pwscf' / 'TiO2_opt' / 'pw.out')


@pytest.fixture(scope='module')
def pwscf_xml_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'pwscf' / 'TiO2_opt' / 'TiO2.save' / 'data-file-schema.xml'
    )


@pytest.fixture(scope='module')
def dos_archive() -> EntryArchive:
    return parse_quantumespresso(DATA_DIR / 'pwscf' / 'W_dos' / 'w.dos.out')


@pytest.fixture(scope='module')
def epw_archive() -> EntryArchive:
    return parse_quantumespresso(DATA_DIR / 'epw' / 'epw.out')


@pytest.fixture(scope='module')
def phonon_archive() -> EntryArchive:
    return parse_quantumespresso(DATA_DIR / 'phonon' / 'ph.out')


@pytest.fixture(scope='module')
def xspectra_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR
        / 'xspectra'
        / 'ms-10734'
        / 'Spectra-1-1-1'
        / '0'
        / 'dipole1'
        / 'xanes.out'
    )


@pytest.fixture(scope='module')
def gipaw_nmr_text_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_out_nmr_out_741' / 'quartz-nmr.out'
    )


@pytest.fixture(scope='module')
def gipaw_nmr_xml_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_xml_nmr_xml' / 'quartz-nmr.xml'
    )


@pytest.fixture(scope='module')
def gipaw_efg_text_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_out_efg_out' / 'quartz-efg.out'
    )


@pytest.fixture(scope='module')
def gipaw_efg_xml_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_xml_efg_xml' / 'quartz-efg.xml'
    )


@pytest.fixture(scope='module')
def gipaw_hyperfine_text_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_out_epr_out' / 'H2O+_hyperfine.out'
    )


@pytest.fixture(scope='module')
def gipaw_hyperfine_xml_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_xml_hyperfine_xml' / 'benzene-hyperfyne.xml'
    )


@pytest.fixture(scope='module')
def gipaw_delta_g_text_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_out_epr_out' / 'H2O+_g-tensor.out'
    )


@pytest.fixture(scope='module')
def gipaw_delta_g_xml_archive() -> EntryArchive:
    return parse_quantumespresso(
        DATA_DIR / 'gipaw' / 'scf_xml_delta_g_xml' / 'benzene-delta_g.xml'
    )
