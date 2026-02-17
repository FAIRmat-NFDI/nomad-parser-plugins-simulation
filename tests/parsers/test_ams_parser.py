from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)

from nomad_simulation_parsers.parsers.ams.parser import AMSParser

LOGGER = get_logger(__name__)
MAINFILE = (
    Path(__file__).resolve().parent.parent
    / 'data'
    / 'ams'
    / 'scf'
    / 'phenylrSmall-metagga.out'
)


def test_parse_file():
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)


def test_workflow_and_scf_steps():
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

    assert isinstance(archive.workflow2, SinglePoint)
    assert isinstance(archive.workflow2.method, SinglePointMethod)
    assert len(archive.workflow2.method.convergence_targets) == 1
    target = archive.workflow2.method.convergence_targets[0]
    assert isinstance(target, EnergyConvergenceTarget)
    assert target.threshold_type == 'absolute'
    assert target.threshold.to('hartree').magnitude == pytest.approx(1e-6)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 1
    assert outputs[0].scf_steps is not None
    assert len(outputs[0].scf_steps.delta_energies_total) == 20
    assert outputs[0].scf_steps.delta_energies_total[0].to(
        'hartree'
    ).magnitude == pytest.approx(0.0)
    assert outputs[0].scf_steps.delta_energies_total[-1].to(
        'hartree'
    ).magnitude == pytest.approx(5.64e-7)
