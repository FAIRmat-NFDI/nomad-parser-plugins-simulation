import numpy as np
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from pytest import approx

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

LOGGER = get_logger(__name__)


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
