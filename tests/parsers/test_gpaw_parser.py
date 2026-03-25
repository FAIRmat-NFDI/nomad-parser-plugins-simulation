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
