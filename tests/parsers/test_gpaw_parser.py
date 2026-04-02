from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)

from nomad_simulation_parsers.parsers.gpaw.parser import GPAWParser

LOGGER = get_logger(__name__)
MAINFILE = Path(__file__).resolve().parent.parent / 'data' / 'gpaw' / 'Fe2.gpw'


def test_parse_file():
    parser = GPAWParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)


def test_workflow_and_scf_steps():
    parser = GPAWParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

    assert isinstance(archive.workflow2, SinglePoint)
    assert isinstance(archive.workflow2.method, SinglePointMethod)
    assert len(archive.workflow2.method.convergence_targets) == 1
    target = archive.workflow2.method.convergence_targets[0]
    assert isinstance(target, EnergyConvergenceTarget)
    assert target.threshold_type == 'absolute'
    assert target.threshold.to('eV').magnitude == pytest.approx(4.0930753554401515e-08)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 1
    assert outputs[0].total_energies is not None
    assert len(outputs[0].total_energies) == 1
    assert outputs[0].total_energies[0].value.to('eV').magnitude == pytest.approx(
        -7.301259879298866
    )
    if outputs[0].electronic_band_structures:
        assert outputs[0].electronic_band_structures[0].value is not None
    if outputs[0].electronic_band_gaps:
        assert outputs[0].electronic_band_gaps[0].value is not None


def test_system_fundamental_quantities_mapping():
    """System gate: parser should populate core model_system quantities used by normalizer."""
    parser = GPAWParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

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
    parser = GPAWParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[0]

    assert output.total_energies or output.total_forces or output.scf_steps is not None

    if output.electronic_band_structures:
        bs = output.electronic_band_structures[0]
        assert bs.value is not None

    if output.electronic_band_gaps:
        assert output.electronic_band_gaps[0].value is not None
