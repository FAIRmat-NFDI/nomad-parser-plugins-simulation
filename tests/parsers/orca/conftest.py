from collections.abc import Generator
from pathlib import Path

import pytest
from nomad import files
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ServerContext
from nomad.utils import create_uuid, get_logger

from nomad_simulation_parsers.parsers.orca.parser import OrcaParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'orca'
LOGGER = get_logger(__name__)


def parse_orca(filename: str) -> EntryArchive:
    archive = EntryArchive()
    OrcaParser().parse(str(DATA_DIR / filename), archive, LOGGER)
    return archive


@pytest.fixture(scope='module')
def ri_mp2_water_archive() -> EntryArchive:
    return parse_orca('RI_MP2_water.out')


@pytest.fixture(scope='module')
def casci_qd_archive() -> EntryArchive:
    return parse_orca('CoPc_CASCI_QD.out')


@pytest.fixture(scope='module')
def dlpno_cc_archive() -> EntryArchive:
    return parse_orca('dlpno-coupled-cluster.out')


class _UploadForTest:
    def __init__(self, upload_id: str, upload_files: files.StagingUploadFiles) -> None:
        self.upload_id = upload_id
        self.upload_files = upload_files


@pytest.fixture(scope='module')
def dft_mos_archive() -> Generator[EntryArchive, None, None]:
    # `dft-print-MOs.out` stores the MO coefficient matrix through the HDF5 backend,
    # which needs an upload-backed `ServerContext`.
    upload_id = f'test_upload_orca_h5_{create_uuid()}'
    upload_files = files.StagingUploadFiles(upload_id, create=True)
    archive = EntryArchive(
        m_context=ServerContext(upload=_UploadForTest(upload_id, upload_files)),
        metadata=EntryMetadata(upload_id=upload_id, entry_id='test_entry_orca_h5'),
    )
    try:
        OrcaParser().parse(str(DATA_DIR / 'dft-print-MOs.out'), archive, LOGGER)
        yield archive
    finally:
        upload_files.delete()
