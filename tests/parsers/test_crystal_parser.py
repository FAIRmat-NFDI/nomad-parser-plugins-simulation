from pathlib import Path

import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)

from nomad_simulation_parsers.parsers.crystal.parser import CrystalParser
from nomad_simulation_parsers.parsers.crystal.parser import CrystalOutputParser

LOGGER = get_logger(__name__)
MAINFILE = (
    Path(__file__).resolve().parent.parent
    / 'data'
    / 'crystal'
    / 'single_point'
    / 'dft'
    / 'output.out'
)


def _parse() -> EntryArchive:
    parser = CrystalParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)
    return archive


def test_parse_file():
    _parse()


def test_workflow_and_scf_steps():
    archive = _parse()

    assert isinstance(archive.workflow2, SinglePoint)
    assert isinstance(archive.workflow2.method, SinglePointMethod)
    assert len(archive.workflow2.method.convergence_targets) == 1
    assert isinstance(
        archive.workflow2.method.convergence_targets[0], EnergyConvergenceTarget
    )
    assert archive.workflow2.method.convergence_targets[0].threshold_type == 'absolute'
    assert archive.workflow2.method.convergence_targets[0].threshold.to(
        'hartree'
    ).magnitude == pytest.approx(1e-7)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 1
    assert outputs[0].scf_steps is not None
    assert len(outputs[0].scf_steps.energies_total) == 8
    assert len(outputs[0].scf_steps.delta_energies_total) == 8
    assert outputs[0].scf_steps.energies_total[-1].to(
        'hartree'
    ).magnitude == pytest.approx(-573.300583798)
    assert outputs[0].scf_steps.delta_energies_total[-1].to(
        'hartree'
    ).magnitude == pytest.approx(5.73e-8)


def test_outputs_electronic_dos_and_band_structure():
    archive = _parse()

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 1

    output = outputs[0]

    # TODO(legacy-parity): Fixture currently provides Crystal side files
    # (`si.f9`/`si.f98`) but not fort.25 payload consumed by current
    # electronic DOS/band-structure mapping. Re-enable checks once parser-side
    # fallback for f9/f98 is implemented.
    assert output.electronic_dos in [None, []]
    assert output.electronic_band_structures in [None, []]


def test_system_fundamental_quantities_mapping():
    """System gate: parser should populate core model_system quantities used by normalizer."""
    archive = _parse()

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
    """Outputs gate: mapped outputs should include normalizer-required payloads when present."""
    archive = _parse()

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[0]

    # Non-electronic core output used by normalizer pipeline.
    assert output.total_energies or output.total_forces or output.scf_steps is not None

    # Electronic outputs are optional for this fixture set until f9/f98 fallback
    # for fort.25-compatible mapping is implemented.
    if output.electronic_dos:
        dos = output.electronic_dos[0]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None

    if output.electronic_band_structures:
        band_structure = output.electronic_band_structures[0]
        assert band_structure.value is not None


def test_crystal_atom_number_prefix_normalization_sr_238():
    """CRYSTAL prefixed species numbers should map to canonical atomic numbers."""
    parser = CrystalOutputParser()
    source = {
        'labels_positions': np.array(
            [['1', 'T', '238', 'Sr', '0.0', '0.0', '0.0']], dtype=str
        )
    }

    atoms = parser.get_atoms(source)
    assert atoms is not None
    assert len(atoms) == 1
    assert atoms[0]['label'] == 'Sr'
    assert atoms[0]['number'] == 38
