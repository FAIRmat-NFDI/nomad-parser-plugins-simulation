from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from pytest import approx

from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)


def test_workflow_convergence_targets():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'GeometryOptimization'
    assert workflow.method is not None
    assert workflow.method.optimization_method is not None
    assert workflow.method.optimization_method.strip() == 'Modified BFGS'

    targets = workflow.method.convergence_targets
    assert targets is not None
    assert len(targets) == 1
    force_target = targets[0]
    assert force_target.m_def.name == 'ForceConvergenceTarget'
    assert force_target.threshold_type == 'maximum'
    assert force_target.threshold.to('eV/angstrom').magnitude == approx(0.01)

    sp_targets = workflow.method.single_point_convergence_targets
    assert sp_targets is not None
    assert len(sp_targets) == 1
    energy_target = sp_targets[0]
    assert energy_target.m_def.name == 'EnergyConvergenceTarget'
    assert energy_target.threshold_type == 'absolute'
    assert energy_target.threshold.to('eV').magnitude == approx(1.0e-6)


def test_scf_steps_quantities():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 5

    expected_n_scf = [11, 7, 8, 6, 7]
    for output, n_scf in zip(outputs, expected_n_scf):
        scf_steps = output.scf_steps
        assert scf_steps is not None
        assert len(scf_steps.delta_energies_total) == n_scf
        assert len(scf_steps.delta_density_rms) == n_scf
        assert len(scf_steps.durations) == n_scf

    # Explicitly check last SCF-step deltas in first geometry-optimization step
    first_steps = outputs[1].scf_steps
    assert first_steps.delta_energies_total[-1].to('eV').magnitude == approx(7.477e-09)
    assert first_steps.delta_density_rms[-1].to('coulomb').magnitude == approx(
        6.375e-08 * 1.602176634e-19
    )


def test_electronic_outputs_mapping():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0

    output = outputs[0]
    if output.electronic_eigenvalues:
        eig = output.electronic_eigenvalues[0]
        assert eig.value is not None
        assert eig.occupation is not None

        # Band-structure compatibility payload should be present when eigenvalues
        # are available from the parser output.
        assert output.electronic_band_structures is not None
        assert len(output.electronic_band_structures) > 0
        sec_bs = output.electronic_band_structures[0]
        assert sec_bs.value is not None

    if output.electronic_dos:
        sec_dos = output.electronic_dos[0]
        assert sec_dos.value is not None
        assert sec_dos.energies is not None
        assert sec_dos.energies.points is not None

    if output.electronic_band_gaps:
        sec_gap = output.electronic_band_gaps[0]
        assert sec_gap.value is not None


def test_system_fundamental_quantities_mapping():
    """System gate for core model_system quantities used by normalizer."""
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    simulation = archive.data
    assert simulation is not None
    assert simulation.model_system is not None
    assert len(simulation.model_system) > 0

    representative = next(
        (s for s in simulation.model_system if getattr(s, 'is_representative', False)),
        simulation.model_system[0],
    )
    assert representative.positions is not None
    assert representative.lattice_vectors is not None
    assert representative.periodic_boundary_conditions is not None

    if representative.particle_states:
        assert all(
            getattr(state, 'chemical_symbol', None) is not None
            for state in representative.particle_states
        )


def test_outputs_contract_for_normalizer():
    """Outputs gate for normalizer-required mapped payloads."""
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[0]

    assert output.total_energies or output.total_forces or output.scf_steps is not None

    if output.electronic_dos:
        dos = output.electronic_dos[0]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None

    if output.electronic_band_gaps:
        assert output.electronic_band_gaps[0].value is not None

    if output.electronic_eigenvalues:
        assert output.electronic_band_structures is not None
        assert output.electronic_band_structures[0].value is not None
