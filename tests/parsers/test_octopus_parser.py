import tempfile
import zipfile
from pathlib import Path

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
import pytest

from nomad_simulation_parsers.parsers.octopus.parser import OctopusParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/stdout.txt', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    assert outputs[-1].total_energies is not None
    assert outputs[-1].total_forces is not None
    assert outputs[-1].electronic_eigenvalues is not None


def test_system_fundamental_quantities_mapping():
    """System gate: parser should populate core model_system quantities used by normalizer."""
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/stdout.txt', archive, LOGGER)

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
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/stdout.txt', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[-1]

    assert output.total_energies or output.total_forces or output.scf_steps is not None

    if output.electronic_eigenvalues:
        eig = output.electronic_eigenvalues[0]
        assert eig.value is not None
        assert eig.occupation is not None

    if output.electronic_band_structures:
        bs = output.electronic_band_structures[0]
        assert bs.value is not None

    if output.electronic_band_gaps:
        assert output.electronic_band_gaps[0].value is not None


def test_root_test_data_octopus_zip_populates_system_from_xyz_sidefile():
    root_dir = Path(__file__).resolve().parents[4]
    zip_path = root_dir / 'test_data' / 'wrZsJFzHT-q4r3MF3H83lA-octopus.zip'
    if not zip_path.is_file():
        pytest.skip(
            'wrZsJFzHT-q4r3MF3H83lA-octopus.zip fixture not available in repository root test_data.'
        )

    parser = OctopusParser()
    archive = EntryArchive()
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmpdir)

        parser.parse(str(Path(tmpdir) / 'output.out'), archive, LOGGER)

    assert archive.data is not None
    assert archive.data.model_system is not None
    assert len(archive.data.model_system) > 0

    representative = next(
        (s for s in archive.data.model_system if getattr(s, 'is_representative', False)),
        archive.data.model_system[0],
    )
    assert representative.positions is not None
    assert representative.lattice_vectors is not None
    assert representative.periodic_boundary_conditions is not None
    assert len(representative.periodic_boundary_conditions) == 3
    assert representative.particle_states is not None
    assert len(representative.particle_states) > 0
