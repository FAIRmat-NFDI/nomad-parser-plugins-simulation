import numpy as np
import pytest

from nomad_simulation_parsers.parsers.fhiaims.common import (
    ControlParser,
    GeometryParser,
)
from nomad_simulation_parsers.parsers.fhiaims.out_parser import FHIAimsOutFileParserLine
from tests.parsers.common import approx


@pytest.mark.unit
class TestGeometryParser:
    def test_reads_geometry_in(self, tmp_path):
        geometry = tmp_path / 'geometry.in'
        geometry.write_text(
            'lattice_vector 2 0 0\n'
            'lattice_vector 0 2 0\n'
            'lattice_vector 0 0 2\n'
            'atom 0 0 0 Si\n'
            'atom_frac 0.5 0.5 0.5 Si\n'
        )

        parser = GeometryParser(str(geometry))
        atoms = parser.get_atoms()

        np.testing.assert_allclose(atoms.cell.array, np.eye(3) * 2)
        assert atoms.get_chemical_symbols() == ['Si', 'Si']
        np.testing.assert_allclose(atoms.get_positions()[1], [1, 1, 1])

    def test_geometry_parser_preserves_fractional_coordinates(self, tmp_path):
        geometry = tmp_path / 'geometry.in'
        geometry.write_text(
            'lattice_vector 4 0 0\n'
            'lattice_vector 0 5 0\n'
            'lattice_vector 0 0 6\n'
            'atom_frac 0.25 0.4 0.5 C\n'
        )

        atoms = GeometryParser(str(geometry)).get_atoms()

        np.testing.assert_allclose(atoms.get_scaled_positions(), [[0.25, 0.4, 0.5]])
        np.testing.assert_allclose(atoms.get_positions(), [[1.0, 2.0, 3.0]])


@pytest.mark.unit
class TestControlParser:
    @pytest.mark.parametrize(
        ('line', 'expected'),
        [('', [0, 0, 0]), ('  k_offset 0.5 0.25 0.0', [0.5, 0.25, 0.0])],
    )
    def test_reads_k_offset_and_applies_gamma_default(self, tmp_path, line, expected):
        control = tmp_path / 'control.in'
        control.write_text('\nk_grid 8 8 8\n' + (line + '\n' if line else ''))
        parser = ControlParser()
        parser.mainfile = str(control)

        assert list(parser.get('k_grid')) == [8, 8, 8]
        offset = parser.get('k_offset')
        if offset is None:
            offset = np.zeros(3)
        np.testing.assert_allclose(offset, expected)

    def test_control_parser_reads_phonon_settings(self, tmp_path):
        control = tmp_path / 'control.in'
        control.write_text(
            '\nphonon supercell 2 3 4\n'
            'phonon displacement 0.01\n'
            'phonon symmetry_thresh 1e-6\n'
        )

        parser = ControlParser()
        parser.mainfile = str(control)

        np.testing.assert_array_equal(parser.get('supercell'), np.diag([2, 3, 4]))
        assert parser.get('displacement') == approx(0.01)
        assert parser.get('symmetry_thresh') == approx(1e-6)


@pytest.mark.unit
class TestFHIAimsOutputReader:
    def test_extracts_structure_from_output(self):
        parser = FHIAimsOutFileParserLine(
            'tests/data/fhiaims/Si_geomopt/out.out'
        )
        parser.line_parsing = True
        parser.allow_overlap = True
        source = parser.to_dict()

        assert {'Number of spin channels': 2} in source['array_size_parameters'][
            'parameter'
        ]
        assert source['structure']['labels'] == ['Si', 'Si']
        np.testing.assert_allclose(
            source['structure']['positions'], [[0, 0, 0], [1.35, 1.35, 1.35]]
        )

    def test_rejects_truncated_input_without_inventing_structure(self, tmp_path):
        output = tmp_path / 'aims.out'
        output.write_text('  Invoking FHI-aims ...\n  Version : 210914\n')
        parser = FHIAimsOutFileParserLine(str(output))
        parser.line_parsing = True
        parser.allow_overlap = True

        assert 'structure' not in parser.to_dict()
