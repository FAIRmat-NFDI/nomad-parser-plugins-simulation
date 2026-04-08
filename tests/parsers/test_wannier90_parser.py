import tempfile
import zipfile
from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.wannier90.parser import (
    Wannier90Parser,
    WInTextParser,
)

LOGGER = get_logger(__name__)


@pytest.fixture
def parsed_archive():
    parser = Wannier90Parser()
    archive = EntryArchive()
    parser.parse('tests/data/wannier90/lco_mlwf/lco.wout', archive, LOGGER)
    return archive


def test_parse_file():
    parser = Wannier90Parser()
    archive = EntryArchive()
    parser.parse('tests/data/wannier90/lco_mlwf/lco.wout', archive, LOGGER)


def test_electronic_outputs_mapping():
    parser = Wannier90Parser()
    archive = EntryArchive()
    parser.parse('tests/data/wannier90/lco_mlwf/lco.wout', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0

    output = outputs[0]
    if output.electronic_dos:
        sec_dos = output.electronic_dos[0]
        assert sec_dos.value is not None
        assert sec_dos.energies is not None
        assert sec_dos.energies.points is not None
        assert sec_dos.energies_origin is not None

    if output.electronic_band_structures:
        sec_band_structure = next(
            (
                sec
                for sec in output.electronic_band_structures
                if getattr(sec, 'value', None) is not None
            ),
            output.electronic_band_structures[0],
        )
        assert sec_band_structure.value is not None
        assert sec_band_structure.highest_occupied is not None

    # Legacy parity: no explicit electronic_band_gaps section for Wannier90.
    assert output.electronic_band_gaps is None or len(output.electronic_band_gaps) == 0


def test_system_fundamental_quantities_mapping(parsed_archive):
    """System gate: parser should populate core model_system quantities used by normalizer."""
    simulation = parsed_archive.data
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


def test_outputs_contract_for_normalizer(parsed_archive):
    """Outputs gate: mapped outputs should include normalizer-required payloads when present."""
    outputs = parsed_archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[0]

    # Wannier90 outputs should provide at least one normalizer-relevant payload.
    assert (
        output.hopping_matrices
        or output.crystal_field_splittings
        or output.electronic_dos
        or output.electronic_band_structures
    )

    if output.electronic_dos:
        dos = output.electronic_dos[0]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None
        assert dos.energies_origin is not None

    if output.electronic_band_structures:
        band_structure = next(
            (
                sec
                for sec in output.electronic_band_structures
                if getattr(sec, 'value', None) is not None
            ),
            output.electronic_band_structures[0],
        )
        assert band_structure.value is not None
        assert band_structure.highest_occupied is not None
        if band_structure.k_path is not None:
            assert getattr(band_structure.k_path, 'points', None) is not None


def test_root_test_data_wannier90_zip_maps_band_structure():
    root_dir = Path(__file__).resolve().parents[4]
    zip_path = root_dir / 'test_data' / '1band-wannier90.zip'
    if not zip_path.is_file():
        pytest.skip('1band-wannier90.zip fixture not available in repository root test_data.')

    parser = Wannier90Parser()
    archive = EntryArchive()
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmpdir)
        parser.parse(str(Path(tmpdir) / '1band.wout'), archive, LOGGER)

    output = archive.data.outputs[0]
    assert output.electronic_band_structures is not None
    assert len(output.electronic_band_structures) > 0
    band_structure = output.electronic_band_structures[0]
    assert band_structure.value is not None
    if band_structure.k_path is not None:
        assert getattr(band_structure.k_path, 'points', None) is not None


class TestWannierQuantumNumberMapping:
    """Test quantum number mapping for Wannier90 orbital symbols."""

    @pytest.fixture
    def parser(self):
        """Create a WInTextParser instance for testing."""
        return WInTextParser()

    def test_s_orbital(self, parser):
        """Test s orbital mapping."""
        ll, ml = parser._get_quantum_numbers_from_symbol('s')
        assert ll == 0
        assert ml == 0

    def test_p_orbitals(self, parser):
        """Test p orbital mappings."""
        assert parser._get_quantum_numbers_from_symbol('px') == (1, -1)
        assert parser._get_quantum_numbers_from_symbol('py') == (1, 0)
        assert parser._get_quantum_numbers_from_symbol('pz') == (1, 1)

    def test_d_orbitals(self, parser):
        """Test d orbital mappings based on Wannier90 Table 3.2."""
        assert parser._get_quantum_numbers_from_symbol('dz2') == (2, 0)
        assert parser._get_quantum_numbers_from_symbol('dxz') == (2, 1)
        assert parser._get_quantum_numbers_from_symbol('dyz') == (2, -1)
        assert parser._get_quantum_numbers_from_symbol('dx2-y2') == (2, 2)
        assert parser._get_quantum_numbers_from_symbol('dxy') == (2, -2)

    def test_f_orbitals(self, parser):
        """Test f orbital mappings."""
        assert parser._get_quantum_numbers_from_symbol('fz3') == (3, 0)
        assert parser._get_quantum_numbers_from_symbol('fxz2') == (3, 1)
        assert parser._get_quantum_numbers_from_symbol('fyz2') == (3, -1)
        assert parser._get_quantum_numbers_from_symbol('fz(x2-y2)') == (3, 2)
        assert parser._get_quantum_numbers_from_symbol('fxyz') == (3, -2)
        assert parser._get_quantum_numbers_from_symbol('fx(x2-3y2)') == (3, 3)
        assert parser._get_quantum_numbers_from_symbol('fy(3x2-y2)') == (3, -3)

    def test_lco_test_case(self, parser):
        """Test the specific orbital used in lco.win test data: Cu:dx2-y2."""
        ll, ml = parser._get_quantum_numbers_from_symbol('dx2-y2')
        assert ll == 2, 'dx2-y2 should have l=2 (d orbital)'
        assert ml == 2, 'dx2-y2 should have ml=+2 according to Wannier90 Table 3.2'

    def test_unknown_symbol(self, parser):
        """Test that unknown symbols return None."""
        result = parser._get_quantum_numbers_from_symbol('unknown')
        assert result is None

    def test_get_orbitals_state_symbol_format(self, parser):
        """Test get_orbitals_state with symbol-based format."""
        states = parser.get_orbitals_state('dx2-y2')
        assert len(states) == 1
        assert states[0]['spin_orbit_state']['l_quantum_number'] == 2
        assert states[0]['spin_orbit_state']['ml_quantum_number'] == 2

    def test_get_orbitals_state_multiple_symbols(self, parser):
        """Test get_orbitals_state with multiple symbols separated by semicolon."""
        states = parser.get_orbitals_state('px;py;pz')
        assert len(states) == 3
        assert states[0]['spin_orbit_state']['l_quantum_number'] == 1
        assert states[0]['spin_orbit_state']['ml_quantum_number'] == -1
        assert states[1]['spin_orbit_state']['ml_quantum_number'] == 0
        assert states[2]['spin_orbit_state']['ml_quantum_number'] == 1
