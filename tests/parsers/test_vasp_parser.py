import numpy as np
import pytest
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ServerContext
from nomad.files import StagingUploadFiles
from nomad.processing import Upload
from nomad.utils import get_logger
from pytest import approx

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

LOGGER = get_logger(__name__)


@pytest.fixture(scope='function')
def test_upload_files():
    return StagingUploadFiles(upload_id='test_upload', create=True)


@pytest.fixture(scope='function')
def test_upload(test_upload_files):
    upload = Upload(upload_id='test_upload')
    # test_upload_files.add_rawfiles('external.h5')
    return upload


@pytest.fixture(scope='function')
def test_context(test_upload):
    return ServerContext(upload=test_upload)


def _parse(mainfile: str) -> EntryArchive:
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse(mainfile, archive, LOGGER)
    return archive


def test_vasprun():
    archive = _parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax')
    assert archive is not None


def test_outcar():
    archive = _parse('tests/data/vasp/AgAc_relax/OUTCAR')
    assert archive is not None


def test_outcar_scf_steps_and_single_point_convergence():
    archive = _parse('tests/data/vasp/AgAc_relax/OUTCAR')
    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'SinglePoint'
    assert workflow.method is not None
    targets = workflow.method.convergence_targets
    assert targets is not None
    assert len(targets) == 1
    assert targets[0].m_def.name == 'EnergyConvergenceTarget'
    assert targets[0].threshold_type == 'absolute'
    assert targets[0].threshold.to('eV').magnitude == approx(1.0e-4)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 1
    scf_steps = outputs[0].scf_steps
    assert scf_steps is not None
    assert len(scf_steps.energies_total) == 11
    assert len(scf_steps.delta_energies_total) == 10
    assert len(scf_steps.durations) == 11
    assert scf_steps.energies_total[-1].to('eV').magnitude == approx(-6.97148118)
    assert np.isfinite(scf_steps.delta_energies_total[-1].to('eV').magnitude)


def test_xml_geometry_optimization_convergence_and_scf_steps():
    archive = _parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax')
    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'GeometryOptimization'
    assert workflow.method is not None

    targets = workflow.method.convergence_targets
    assert targets is not None
    assert len(targets) == 1
    assert targets[0].m_def.name == 'EnergyConvergenceTarget'
    assert targets[0].threshold_type == 'absolute'
    assert targets[0].threshold.to('eV').magnitude == approx(1.0e-3)

    sp_targets = workflow.method.single_point_convergence_targets
    assert sp_targets is not None
    assert len(sp_targets) == 1
    assert sp_targets[0].m_def.name == 'EnergyConvergenceTarget'
    assert sp_targets[0].threshold_type == 'absolute'
    assert sp_targets[0].threshold.to('eV').magnitude == approx(1.0e-4)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 3
    expected_scf_counts = [12, 10, 6]
    for output, n_scf in zip(outputs, expected_scf_counts):
        scf_steps = output.scf_steps
        assert scf_steps is not None
        assert len(scf_steps.energies_total) == n_scf
        assert len(scf_steps.delta_energies_total) == n_scf - 1
        assert len(scf_steps.durations) == n_scf


def test_chgcar(test_context, test_upload):
    archive = EntryArchive(
        m_context=test_context,
        metadata=EntryMetadata(upload_id=test_upload.upload_id, entry_id='test_entry'),
    )
    parser = VASPParser()
    parser.parse('tests/data/vasp/with_chgcar/OUTCAR', archive, LOGGER)
    assert len(archive.data.outputs) == 14
    outputs = archive.data.outputs[-1]
    assert len(outputs.charge_density) == 1
    assert outputs.charge_density[0].value_h5_dataset is not None
    with outputs.charge_density[0].value_h5_dataset as dataset:
        dataset_array = dataset[:]
        assert np.shape(dataset_array) == (10, 10, 10)
        assert dataset_array[9][9][9] == approx(0.18013097030e-05)
