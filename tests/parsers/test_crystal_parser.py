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

    assert len(output.electronic_dos or []) > 0
    sec_dos = output.electronic_dos[0]
    assert sec_dos.value is not None
    assert sec_dos.energies is not None
    assert sec_dos.energies.points is not None

    assert len(output.electronic_band_structures or []) > 0
    sec_band_structure = output.electronic_band_structures[0]
    assert sec_band_structure.value is not None
