import os
from pathlib import Path

import numpy as np
import pytest

from nomad_simulation_parsers.parsers.gromacs.edr_parser import GromacsEDRParser
from nomad_simulation_parsers.parsers.gromacs.log_parser import GromacsLogParser
from nomad_simulation_parsers.parsers.gromacs.mdanalysis_parser import (
    GromacsMDAnalysisParser as GromacsMDAnalysisFileParser,
)
from nomad_simulation_parsers.parsers.gromacs.mdanalysis_parser import (
    group_indices_by_value,
)
from nomad_simulation_parsers.parsers.gromacs.mdp_parser import GromacsMdpParser
from nomad_simulation_parsers.parsers.gromacs.xvg_parser import GromacsXvgParser
from tests.parsers.common import assert_approx

LARGE_DATA_DIR = (
    Path(os.environ['NOMAD_SIM_PARSERS_LARGE_FIXTURE_ROOT']) / 'gromacs'
    if 'NOMAD_SIM_PARSERS_LARGE_FIXTURE_ROOT' in os.environ
    else Path(__file__).resolve().parents[2] / 'data' / 'gromacs'
)


def system_node(name, branch_label, particle_indices, sub_systems=None):
    node = {
        'name': name,
        'branch_label': branch_label,
        'particle_indices': np.array(particle_indices),
    }
    if sub_systems:
        node['sub_systems'] = sub_systems
    return node


def test_group_indices_by_value_groups_only_selected_indices_once():
    values = np.array(['B', 'A', 'B', 'A', 'B', 'A'])

    grouped = group_indices_by_value([5, 2, 3, 0], values)

    assert list(grouped) == ['A', 'B']
    assert grouped['A'].tolist() == [3, 5]
    assert grouped['B'].tolist() == [0, 2]


@pytest.mark.unit
class TestSystemNameDisambiguation:
    def test_conflicting_names_are_disambiguated_and_retained(self):
        repeated_large = system_node('PEO', 'molecule', range(8))
        same_signature = system_node('PEO', 'molecule', range(8, 16))
        repeated_small = system_node('PEO', 'monomer', range(3))
        root = system_node(
            'group_PEO',
            'molecule_group',
            range(16),
            [repeated_large, same_signature, repeated_small],
        )

        GromacsMDAnalysisFileParser.disambiguate_system_names([root])

        assert root['name'] == 'group_PEO'
        assert repeated_small['name'] == 'PEO_1'
        assert repeated_large['name'] == 'PEO'
        assert same_signature['name'] == repeated_large['name']
        assert root['sub_systems'] == [
            repeated_large,
            same_signature,
            repeated_small,
        ]
        assert repeated_small['particle_indices'].tolist() == [0, 1, 2]
        assert repeated_large['particle_indices'].tolist() == list(range(8))

    def test_same_type_and_size_names_are_left_unchanged(self):
        first = system_node('PEO', 'molecule', range(3))
        second = system_node('PEO', 'molecule', range(3, 6))

        GromacsMDAnalysisFileParser.disambiguate_system_names([first, second])

        assert first['name'] == second['name'] == 'PEO'

    def test_same_sized_groups_under_different_chains_get_distinct_names(self):
        chains = []
        for chain_index, chain_name in enumerate(
            ['Protein_chain_A', 'Protein_chain_B']
        ):
            chain_start = chain_index * 6
            residues = [
                system_node('SA', 'monomer', range(start, start + 3))
                for start in (chain_start, chain_start + 3)
            ]
            monomer_group = system_node(
                'group_SA',
                'monomer_group',
                range(chain_start, chain_start + 6),
                residues,
            )
            molecule = system_node(
                chain_name,
                'molecule',
                range(chain_start, chain_start + 6),
                [monomer_group],
            )
            chains.append(
                system_node(
                    f'group_{chain_name}',
                    'molecule_group',
                    range(chain_start, chain_start + 6),
                    [molecule],
                )
            )

        GromacsMDAnalysisFileParser.disambiguate_system_names(chains)

        chain_a_group = chains[0]['sub_systems'][0]['sub_systems'][0]
        chain_b_group = chains[1]['sub_systems'][0]['sub_systems'][0]
        assert chain_a_group['name'] == 'group_SA'
        assert chain_b_group['name'] == 'group_SA_1'
        assert [group['name'] for group in chain_a_group['sub_systems']] == [
            'SA',
            'SA',
        ]
        assert [group['name'] for group in chain_b_group['sub_systems']] == [
            'SA_1',
            'SA_1',
        ]
        assert chain_a_group['particle_indices'].tolist() == list(range(6))
        assert chain_b_group['particle_indices'].tolist() == list(range(6, 12))

    def test_disambiguated_names_avoid_existing_names(self):
        first = system_node('PEO', 'molecule', range(2))
        other_name = system_node('PEO_1', 'molecule', range(2, 4))
        conflicting = system_node('PEO', 'monomer', [0])

        GromacsMDAnalysisFileParser.disambiguate_system_names(
            [first, other_name, conflicting]
        )

        assert first['name'] == 'PEO'
        assert other_name['name'] == 'PEO_1'
        assert conflicting['name'] == 'PEO_2'


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

    @pytest.mark.large_fixture
    @pytest.mark.parametrize(
        ('directory', 'topology', 'trajectory', 'n_atoms'),
        [
            ('protein_fsfg', 'nvt.tpr', 'nvt.trr', 45194),
            ('cgwater', 'cgwater.tpr', 'cgwater.trr', 1000),
        ],
    )
    def test_reads_large_trajectory_sets(
        self, directory, topology, trajectory, n_atoms
    ):
        base = LARGE_DATA_DIR / directory
        parser = GromacsMDAnalysisFileParser()
        parser.mainfile = str(base / topology)
        parser.auxilliary_files = [str(base / trajectory)]

        assert parser.get_n_atoms(0) == n_atoms
        assert parser.get_frame_data(0)['positions'].shape == (n_atoms, 3)


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
