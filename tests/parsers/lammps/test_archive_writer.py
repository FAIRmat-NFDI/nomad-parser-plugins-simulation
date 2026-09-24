import pytest

from nomad_simulation_parsers.parsers.lammps.parser import LammpsArchiveWriter
from tests.parsers.common import approx


@pytest.mark.unit
class TestLammpsArchiveWriter:
    def test_archive_writer_converts_values_to_log_units(self):
        writer = LammpsArchiveWriter()
        writer._log_parser._units = {'distance': 1, 'energy': 1}

        assert writer.apply_unit(2, 'distance') == approx(2)
        assert writer.apply_unit(3, 'energy') == approx(3)
