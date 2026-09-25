from pathlib import Path

import pytest

from nomad_simulation_parsers.parsers.gromacs.edr_parser import GromacsEDRParser
from nomad_simulation_parsers.parsers.gromacs.log_parser import GromacsLogParser
from nomad_simulation_parsers.parsers.gromacs.mdanalysis_parser import (
    GromacsMDAnalysisParser as GromacsMDAnalysisFileParser,
)
from nomad_simulation_parsers.parsers.gromacs.mdp_parser import GromacsMdpParser
from nomad_simulation_parsers.parsers.gromacs.xvg_parser import GromacsXvgParser
from tests.parsers.common import assert_approx


@pytest.mark.unit
class TestGromacsLogReader:
    @pytest.fixture
    def parser(self):
        return GromacsLogParser()

    def test_reads_input_parameter_quantities(self, tmp_path, parser):
        log = tmp_path / 'md.log'
        log.write_text(
            'GROMACS version: 2025.0\n'
            'Input Parameters:\n'
            ' integrator: md\n'
            ' pbc: xyz\n'
            ' constraints: all-bonds\n'
            ' ref-t[0] = {300,300}\n\n'
        )
        parser.mainfile = str(log)
        parser.parse()

        assert parser.results['version'] == 2025.0
        parameters = parser.results['input_parameters']
        assert parameters['integrator'] == 'md'
        assert parameters['pbc'] == 'xyz'
        assert parameters['constraints'] == 'all-bonds'
        assert_approx(parameters['ref-t'], [[300, 300]])


@pytest.mark.unit
class TestGromacsMDPReader:
    @pytest.fixture
    def parser(self):
        return GromacsMdpParser()

    def test_reads_typed_mdp_quantities(self, tmp_path, parser):
        mdp = tmp_path / 'md.mdp'
        mdp.write_text('integrator = md\nnsteps = 1000\nref-t[0] = 300, 300\n')
        parser.mainfile = str(mdp)
        parser.parse()

        parameters = parser.results['input_parameters']
        assert parameters['integrator'] == 'md'
        assert parameters['nsteps'] == 1000
        assert_approx(parameters['ref-t'], [[300, 300]])


@pytest.mark.unit
class TestGromacsXVGReader:
    @pytest.fixture
    def parser(self):
        return GromacsXvgParser()

    def test_reads_metadata_and_values(self, tmp_path, parser):
        xvg = tmp_path / 'energy.xvg'
        xvg.write_text(
            '@    title "Potential"\n'
            '@    xaxis  label "Time (ps)"\n'
            '@    yaxis  label "Energy (kJ/mol)"\n'
            '@ s0 legend "Potential"\n'
            '0.0 -10.0\n'
            '1.0 -9.5\n'
        )
        parser.mainfile = str(xvg)
        parser.parse()

        assert parser.results['title'] == 'Potential'
        assert parser.results['column_headers'] == ['Potential']
        assert_approx(parser.results['column_vals'], [[0.0, -10.0], [1.0, -9.5]])


@pytest.mark.unit
class TestGromacsMDAnalysisReader:
    @pytest.fixture
    def parser(self):
        base = Path(__file__).resolve().parents[2] / 'data' / 'gromacs' / 'water'
        parser = GromacsMDAnalysisFileParser()
        parser.mainfile = str(base / 'reference_s.tpr')
        parser.auxilliary_files = [str(base / 'reference_s.trr')]
        return parser

    def test_reads_topology_and_trajectory_quantities(self, parser):
        assert parser.get_n_atoms(0) == 648
        assert parser.get_atom_labels(0)[:3] == ['O', 'H', 'H']

        frame = parser.get_frame_data(0)
        assert frame['positions'].shape == (648, 3)
        assert frame['velocities'].shape == (648, 3)
        assert frame['lattice_vectors'].shape == (3, 3)


@pytest.mark.unit
class TestGromacsEDRReader:
    @pytest.fixture
    def parser(self):
        base = Path(__file__).resolve().parents[2] / 'data' / 'gromacs' / 'water'
        parser = GromacsEDRParser()
        parser.mainfile = str(base / 'reference_s.edr')
        return parser

    def test_reads_energy_quantities(self, parser):
        assert 'Temperature' in parser.keys()
        assert 'Potential' in parser.keys()

        parser.parse('Temperature')

        assert len(parser.results['Temperature']) == parser.length
        assert_approx(
            parser.results['Temperature'][:3],
            [11.41081619, 509.01992798, 507.32138062],
            atol=1e-6,
        )
