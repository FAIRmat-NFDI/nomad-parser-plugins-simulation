import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.wannier90.parser import (
    Wannier90Parser,
    WInTextParser,
)

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = Wannier90Parser()
    archive = EntryArchive()
    parser.parse('tests/data/wannier90/lco_mlwf/lco.wout', archive, LOGGER)


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
