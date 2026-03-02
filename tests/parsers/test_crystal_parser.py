from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)

from nomad_simulation_parsers.parsers.crystal.parser import CrystalParser

LOGGER = get_logger(__name__)
MAINFILE = (
    Path(__file__).resolve().parent.parent
    / 'data'
    / 'crystal'
    / 'single_point'
    / 'dft'
    / 'output.out'
)


def test_parse_file():
    parser = CrystalParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)


def test_workflow_and_scf_steps():
    parser = CrystalParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

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
