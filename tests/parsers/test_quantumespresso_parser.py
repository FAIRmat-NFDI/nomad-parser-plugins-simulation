import tempfile
import zipfile
from pathlib import Path

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from pytest import approx
import pytest

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoParser,
)

LOGGER = get_logger(__name__)


def _parse(mainfile: str) -> EntryArchive:
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(mainfile, archive, LOGGER)
    return archive


def test_pwscf():
    archive = _parse('tests/data/quantumespresso/pwscf/TiO2_opt/pw.out')
    assert archive is not None


def test_pwscf_xml():
    archive = _parse(
        'tests/data/quantumespresso/pwscf/TiO2_opt/TiO2.save/data-file-schema.xml',
    )
    assert archive is not None


def test_epw():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/epw/epw.out', archive, LOGGER)


def test_phonon():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/phonon/ph.out', archive, LOGGER)


def test_xspectra():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/xspectra/ms-10734/Spectra-1-1-1/0/dipole1/xanes.out',
        archive,
        LOGGER,
    )


def test_gipaw_nmr_text():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_out_nmr_out_741/quartz-nmr.out',
        archive,
        LOGGER,
    )


def test_gipaw_nmr_xml():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_xml_nmr_xml/quartz-nmr.xml',
        archive,
        LOGGER,
    )


def test_gipaw_efg_text():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_out_efg_out/quartz-efg.out',
        archive,
        LOGGER,
    )


def test_gipaw_efg_xml():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_xml_efg_xml/quartz-efg.xml',
        archive,
        LOGGER,
    )


def test_gipaw_epr_hyperfine_text():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_out_epr_out/H2O+_hyperfine.out',
        archive,
        LOGGER,
    )


def test_gipaw_epr_hyperfine_xml():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_xml_hyperfine_xml/benzene-hyperfyne.xml',
        archive,
        LOGGER,
    )


def test_gipaw_epr_deltag_text():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_out_epr_out/H2O+_g-tensor.out',
        archive,
        LOGGER,
    )


def test_gipaw_epr_deltag_xml():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_xml_delta_g_xml/benzene-delta_g.xml',
        archive,
        LOGGER,
    )


def test_pwscf_workflow_and_scf_steps():
    archive = _parse('tests/data/quantumespresso/pwscf/TiO2_opt/pw.out')
    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'GeometryOptimization'
    assert workflow.method is not None

    # In this fixture we have reliable SCF threshold but no parsed ionic thresholds.
    assert workflow.method.convergence_targets in [None, []]
    sp_targets = workflow.method.single_point_convergence_targets
    assert sp_targets is not None
    assert len(sp_targets) == 1
    assert sp_targets[0].m_def.name == 'EnergyConvergenceTarget'
    assert sp_targets[0].threshold_type == 'absolute'
    assert sp_targets[0].threshold.to('rydberg').magnitude == approx(1.0e-8)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 2
    assert outputs[0].scf_steps is not None
    assert outputs[1].scf_steps is not None
    assert len(outputs[0].scf_steps.energies_total) == 12
    assert len(outputs[1].scf_steps.energies_total) == 14
    assert len(outputs[0].scf_steps.delta_energies_total) == 12
    assert len(outputs[1].scf_steps.delta_energies_total) == 14
    assert outputs[0].electronic_eigenvalues is not None
    assert outputs[0].electronic_eigenvalues[0].occupation is not None


def test_pwscf_xml_workflow_and_scf_steps():
    archive = _parse(
        'tests/data/quantumespresso/pwscf/TiO2_opt/TiO2.save/data-file-schema.xml'
    )
    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'GeometryOptimization'
    assert workflow.method is not None

    targets = workflow.method.convergence_targets
    assert targets is not None
    assert len(targets) == 2
    targets_by_name = {target.m_def.name: target for target in targets}
    assert targets_by_name['ForceConvergenceTarget'].threshold_type == 'maximum'
    assert targets_by_name['ForceConvergenceTarget'].threshold.to(
        'rydberg/bohr'
    ).magnitude == approx(2.5e-4)
    assert targets_by_name['EnergyConvergenceTarget'].threshold_type == 'absolute'
    assert targets_by_name['EnergyConvergenceTarget'].threshold.to(
        'rydberg'
    ).magnitude == approx(5.0e-5)

    sp_targets = workflow.method.single_point_convergence_targets
    assert sp_targets is not None
    assert len(sp_targets) == 1
    assert sp_targets[0].threshold.to('rydberg').magnitude == approx(5.0e-9)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 2
    for output in outputs:
        assert output.scf_steps is not None
        assert len(output.scf_steps.energies_total) == 5
        assert len(output.scf_steps.delta_energies_total) == 5


def test_root_test_data_pwscf_dos_zip_populates_system_and_dos():
    root_dir = Path(__file__).resolve().parents[4]
    zip_path = root_dir / 'test_data' / 'DOS-quantumespresso.zip'
    if not zip_path.is_file():
        pytest.skip('DOS-quantumespresso.zip fixture not available in repository root test_data.')

    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmpdir)

        archive = _parse(str(Path(tmpdir) / 'W.out'))

    simulation = archive.data
    assert simulation is not None
    assert simulation.model_system is not None
    assert len(simulation.model_system) > 0

    representative_candidates = [
        s for s in simulation.model_system if s.is_representative
    ]
    representative = (
        representative_candidates[0]
        if representative_candidates
        else simulation.model_system[0]
    )
    assert representative.positions is not None
    assert representative.lattice_vectors is not None
    assert representative.periodic_boundary_conditions is not None
    assert len(representative.periodic_boundary_conditions) == 3

    outputs = simulation.outputs
    assert outputs is not None
    assert len(outputs) > 0
    for output in outputs:
        assert output.model_system_ref is not None

    output = outputs[0]
    assert output.electronic_dos is not None
    assert len(output.electronic_dos) > 0
    dos = output.electronic_dos[0]
    assert dos.value is not None
    assert dos.energies is not None
    assert dos.energies.points is not None
