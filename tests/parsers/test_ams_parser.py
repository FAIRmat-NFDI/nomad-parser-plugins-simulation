import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)

from nomad_simulation_parsers.parsers.ams.parser import AMSParser, MainfileParser

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


# TODO: Remove this skip once EntryMetadata.auxiliary_files is available.
@pytest.mark.skip(reason='requires EntryMetadata.auxiliary_files in nomad.datamodel')
def test_parse_file_records_parsed_blocks():
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

    assert archive.metadata.auxiliary_files
    assert any(file.parsed_blocks for file in archive.metadata.auxiliary_files)


def test_periodic_boundary_conditions_mapping_input_variants():
    parser = MainfileParser()

    assert parser.get_periodic_boundary_conditions({'lattice_vectors': None}) == [
        False,
        False,
        False,
    ]
    assert parser.get_periodic_boundary_conditions(np.eye(3) * ureg.angstrom) == [
        True,
        True,
        True,
    ]
    assert parser.get_periodic_boundary_conditions(np.eye(2, 3) * ureg.angstrom) == [
        True,
        True,
        False,
    ]


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

    assert outputs[0].electronic_band_gaps is not None
    assert len(outputs[0].electronic_band_gaps) == 1
    gap = outputs[0].electronic_band_gaps[0].value.to('hartree').magnitude
    assert gap == pytest.approx(0.097, rel=1e-3)


def test_system_fundamental_quantities_mapping():
    """System gate for core model_system quantities used by normalizer."""
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

    simulation = archive.data
    assert simulation is not None
    assert simulation.model_system is not None
    assert len(simulation.model_system) > 0

    representative = next(
        (s for s in simulation.model_system if s.is_representative),
        simulation.model_system[0],
    )
    assert representative.positions is not None
    if representative.periodic_boundary_conditions is not None:
        assert len(representative.periodic_boundary_conditions) == 3
        assert all(
            isinstance(flag, bool)
            for flag in representative.periodic_boundary_conditions
        )

    # This fixture is molecular, so lattice vectors can be absent.
    # If lattice vectors are present, ensure periodic axes are available.
    if representative.lattice_vectors is not None:
        assert representative.periodic_boundary_conditions is not None

    if representative.particle_states:
        assert all(
            state.chemical_symbol is not None
            for state in representative.particle_states
        )


def test_outputs_contract_for_normalizer():
    """Outputs gate for normalizer-required mapped payloads."""
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse(str(MAINFILE), archive, LOGGER)

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


def test_root_test_data_ams_zip_outputs_and_system_links():
    root_dir = Path(__file__).resolve().parents[4]
    zip_path = root_dir / 'test_data' / 'FAJZIC_fair_op-ams.zip'
    if not zip_path.is_file():
        pytest.skip(
            'FAJZIC_fair_op-ams.zip fixture not available in repository root test_data.'
        )

    parser = AMSParser()
    archive = EntryArchive()
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmpdir)

        parser.parse(str(Path(tmpdir) / 'FAJZIC_fair_op.out'), archive, LOGGER)

    simulation = archive.data
    assert simulation is not None
    assert simulation.model_system is not None
    assert len(simulation.model_system) > 0
    assert simulation.outputs is not None
    assert len(simulation.outputs) > 0

    assert all(out.model_system_ref is not None for out in simulation.outputs)
    step_outputs = [
        out for out in simulation.outputs if 'step' in out.m_def.all_quantities
    ]
    assert all(out.step is not None for out in step_outputs)

    dos_outputs = [out for out in simulation.outputs if out.electronic_dos]
    assert len(dos_outputs) > 0
    dos = dos_outputs[0].electronic_dos[0]
    assert dos.value is not None
    assert dos.energies is not None
    assert dos.energies.points is not None


@pytest.mark.parametrize(
    'energies_channels, occupations_channels, expected',
    [
        pytest.param(
            [[-5.0, -3.0, 1.0, 2.0]],
            [[[1.0, 1.0, 0.0, 0.0]]],
            [dict(value=4.0, spin_channel=0)],
            id='single-spin-semiconductor',
        ),
        pytest.param(
            [[-5.0, -3.0, 1.0, 2.0], [-4.0, -2.0, 0.5, 3.0]],
            [[[1.0, 1.0, 0.0, 0.0]], [[1.0, 1.0, 0.0, 0.0]]],
            [dict(value=4.0, spin_channel=0), dict(value=2.5, spin_channel=1)],
            id='two-spin-channels',
        ),
        pytest.param(
            [[-5.0, -3.0, 1.0]],
            [[[1.0, 1.0]]],
            [],
            id='shape-mismatch-skipped',
        ),
        pytest.param(
            [[-5.0, -3.0, 1.0]],
            [[]],
            [],
            id='empty-occupations-skipped',
        ),
    ],
)
def test_band_gaps_from_tuple_payload(
    energies_channels, occupations_channels, expected
):
    """Band gaps from the tuple/list fallback payload of band-energy ranges.

    Regression test: this path previously called the band-gap utility with an
    unknown keyword and indexed its dict return as a list, raising at runtime.
    """
    parser = MainfileParser()

    gaps = parser.get_band_gaps((energies_channels, None, occupations_channels))

    assert len(gaps) == len(expected)
    for gap, ref in zip(gaps, expected):
        assert gap['value'] == pytest.approx(ref['value'])
        assert gap['spin_channel'] == ref['spin_channel']
