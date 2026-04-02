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


def test_vasprun_system_and_electronic_outputs():
    archive = _parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax')

    sec_system = archive.data.model_system[0]
    assert sec_system.positions is not None
    assert sec_system.lattice_vectors is not None
    assert sec_system.periodic_boundary_conditions == [True, True, True]
    assert len(sec_system.particle_states) == len(sec_system.positions)
    assert sec_system.particle_states[0].chemical_symbol == 'Ac'

    sec_output = archive.data.outputs[0]
    assert len(sec_output.electronic_band_structures) > 0
    assert len(sec_output.electronic_band_gaps) > 0
    assert len(sec_output.electronic_dos) > 0

    sec_band_structure = sec_output.electronic_band_structures[0]
    assert sec_band_structure.value is not None

    sec_band_gap = sec_output.electronic_band_gaps[0]
    assert sec_band_gap.value is not None

    sec_dos = sec_output.electronic_dos[0]
    assert sec_dos.value is not None
    assert sec_dos.energies is not None
    assert sec_dos.energies.points is not None


def test_outcar():
    archive = _parse('tests/data/vasp/AgAc_relax/OUTCAR')
    assert archive is not None

    sec_system = archive.data.model_system[0]
    assert sec_system.positions is not None
    assert sec_system.lattice_vectors is not None
    assert sec_system.periodic_boundary_conditions == [True, True, True]


def test_outcar_electronic_outputs_from_doscar_and_eigenvalues():
    archive = _parse('tests/data/vasp/AgAc_relax/OUTCAR')

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[0]

    assert output.electronic_band_structures is not None
    assert len(output.electronic_band_structures) > 0
    band_structure = output.electronic_band_structures[0]
    assert band_structure.value is not None

    if output.electronic_band_gaps:
        assert output.electronic_band_gaps[0].value is not None

    assert output.electronic_dos is not None
    assert len(output.electronic_dos) > 0
    dos = output.electronic_dos[0]
    assert dos.value is not None
    assert dos.energies is not None
    assert dos.energies.points is not None


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
