from types import SimpleNamespace

import numpy as np
import pytest
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem

from nomad_simulation_parsers.parsers.gromacs.mdanalysis_parser import (
    GromacsMDAnalysisParser,
)
from nomad_simulation_parsers.parsers.lammps.parser import LammpsArchiveWriter
from nomad_simulation_parsers.parsers.utils.molecular_hierarchy import (
    build_molecular_hierarchy,
    disambiguate_system_names,
    group_indices_by_value,
    infer_residue_types,
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


def interleaved_polymer_arrays():
    """Two eight-site chains, each with two end groups and two interior residues."""
    types = np.empty(16, dtype=object)
    bonds_by_particle = [[] for _ in range(16)]
    for molecule in (0, 1):
        indices = list(range(molecule, 16, 2))
        types[indices] = list('ABBBBBBC')
        for first, second in zip(indices[:-1], indices[1:], strict=True):
            bonds_by_particle[first].append((second, '1'))
            bonds_by_particle[second].append((first, '1'))
    arrays = {
        'moltypes': np.array(['P'] * 16, dtype=object),
        'molnums': np.array([0, 1] * 8),
        'resids': np.repeat([1, 2, 3, 4], 4),
        'resnames': np.repeat([1, 2, 3, 4], 4).astype(object),
        'types': types,
        'elements': np.array(['C'] * 16),
        'labels': np.array(['C'] * 16),
    }
    return arrays, bonds_by_particle


def test_unlabeled_interleaved_residues_share_types_not_occurrence_names():
    arrays, bonds = interleaved_polymer_arrays()
    arrays['resnames'] = infer_residue_types(arrays, bonds)
    group = build_molecular_hierarchy(arrays)['P']

    molecules = group['sub_systems']
    assert len(molecules) == 2
    assert [molecule['name'] for molecule in molecules] == ['P', 'P']
    for molecule_index, molecule in enumerate(molecules):
        np.testing.assert_array_equal(
            molecule['particle_indices'], np.arange(molecule_index, 16, 2)
        )
        groups = molecule['sub_systems']
        assert [subgroup['name'] for subgroup in groups] == [
            'group_1',
            'group_2',
            'group_3',
        ]
        assert [
            [monomer['name'] for monomer in subgroup['sub_systems']]
            for subgroup in groups
        ] == [['1'], ['2', '2'], ['3']]
        assert [len(subgroup['particle_indices']) for subgroup in groups] == [2, 4, 2]
        assert all(
            set(subgroup['particle_indices']) <= set(molecule['particle_indices'])
            for subgroup in groups
        )


def test_residue_type_inference_keeps_different_bond_topologies_separate():
    bonds = [[] for _ in range(6)]
    for first, second in [(0, 1), (1, 2), (3, 4), (3, 5)]:
        bonds[first].append((second, '1'))
        bonds[second].append((first, '1'))
    arrays = {
        'moltypes': np.array(['P'] * 6),
        'molnums': np.zeros(6, dtype=int),
        'resids': np.repeat([1, 2], 3),
        'resnames': np.repeat([1, 2], 3),
        'types': np.array(['A', 'B', 'C'] * 2),
    }

    names = infer_residue_types(arrays, bonds)

    assert names.tolist() == ['1'] * 3 + ['2'] * 3


def test_missing_typed_bonds_preserve_residue_labels():
    arrays, _ = interleaved_polymer_arrays()
    original = arrays['resnames'].copy()

    names = infer_residue_types(arrays, None)

    np.testing.assert_array_equal(names, original)
    np.testing.assert_array_equal(arrays['resnames'], original)


def test_builder_disambiguates_types_without_renaming_occurrences():
    arrays, bonds = interleaved_polymer_arrays()
    arrays['moltypes'][:] = '1'
    arrays['resnames'] = infer_residue_types(arrays, bonds)
    group = build_molecular_hierarchy(arrays)['1']

    assert group['name'] == 'group_1_0'
    molecules = group['sub_systems']
    assert [molecule['name'] for molecule in molecules] == ['1_0', '1_0']
    for molecule in molecules:
        assert molecule['sub_systems'][0]['name'] == 'group_1_1'
        assert molecule['sub_systems'][0]['sub_systems'][0]['name'] == '1_1'


def test_composition_uses_final_child_names_and_preserves_leaf_formulas():
    arrays, bonds = interleaved_polymer_arrays()
    arrays['moltypes'][:] = 'PEO'
    residue_types = infer_residue_types(arrays, bonds)
    arrays['resnames'] = np.array(
        ['PEO' if name == '2' else 'CH3' for name in residue_types], dtype=object
    )

    group = build_molecular_hierarchy(
        arrays, chemical_formula=lambda indices: 'CustomLeaf(1)'
    )['PEO']

    assert group['name'] == 'group_PEO_0'
    assert group['composition_formula'] == 'PEO_0(2)'
    for molecule in group['sub_systems']:
        assert molecule['name'] == 'PEO_0'
        assert molecule['composition_formula'] == 'group_CH3(1)group_PEO_1(1)'
        by_name = {child['name']: child for child in molecule['sub_systems']}
        assert by_name['group_CH3']['composition_formula'] == 'CH3(2)'
        assert by_name['group_PEO_1']['composition_formula'] == 'PEO_1(2)'
        for subgroup in molecule['sub_systems']:
            assert all(
                monomer['composition_formula'] == 'CustomLeaf(1)'
                for monomer in subgroup['sub_systems']
            )


def test_builder_rejects_misaligned_particle_arrays():
    arrays, _ = interleaved_polymer_arrays()
    arrays['resids'] = arrays['resids'][:-1]

    with pytest.raises(ValueError, match='equal particle counts'):
        build_molecular_hierarchy(arrays)


def test_gromacs_and_lammps_build_the_same_named_hierarchy(monkeypatch):
    arrays, bonds = interleaved_polymer_arrays()
    arrays['resnames'] = infer_residue_types(arrays, bonds)
    reader = GromacsMDAnalysisParser()
    monkeypatch.setattr(
        reader,
        'get',
        lambda key, default=None: arrays if key == 'atoms_info' else default,
    )
    gromacs = reader.parse_molecular_hierarchy()['P']
    writer = LammpsArchiveWriter()
    monkeypatch.setattr(
        writer._mdanalysistraj_parser,
        'get',
        lambda key, default=None: arrays if key == 'atoms_info' else default,
    )
    simulation = Simulation(model_system=[ModelSystem()])
    writer._parse_molecular_hierarchy(simulation)
    lammps = simulation.model_system[0].sub_systems[0]

    def compare(gromacs_node, lammps_node):
        assert gromacs_node['name'] == lammps_node.name
        assert gromacs_node['branch_label'] == lammps_node.branch_label
        np.testing.assert_array_equal(
            gromacs_node['particle_indices'], lammps_node.particle_indices
        )
        children = gromacs_node.get('sub_systems', [])
        assert len(children) == len(lammps_node.sub_systems)
        for gromacs_child, lammps_child in zip(
            children, lammps_node.sub_systems, strict=True
        ):
            compare(gromacs_child, lammps_child)

    compare(gromacs, lammps)


@pytest.mark.parametrize('typed_bonds', [True, False])
def test_lammps_fallback_names_identify_roles_and_types(monkeypatch, typed_bonds):
    arrays, bonds = interleaved_polymer_arrays()
    writer = LammpsArchiveWriter()
    writer._mdanalysistraj_parser = SimpleNamespace(
        universe=SimpleNamespace(atoms=SimpleNamespace()),
        get=lambda key, default=None: arrays if key == 'atoms_info' else default,
    )
    monkeypatch.setattr(
        writer,
        '_typed_bonds_by_particle',
        lambda universe, types: bonds if typed_bonds else None,
    )
    simulation = Simulation(model_system=[ModelSystem()])

    writer._parse_molecular_hierarchy(simulation)

    group = simulation.model_system[0].sub_systems[0]
    assert group.name == 'group_molecule_1'
    assert [molecule.name for molecule in group.sub_systems] == [
        'molecule_1',
        'molecule_1',
    ]
    expected_names = (
        ['monomer_1', 'monomer_2', 'monomer_3']
        if typed_bonds
        else ['monomer_1', 'monomer_2', 'monomer_3', 'monomer_4']
    )
    expected_counts = [1, 2, 1] if typed_bonds else [1, 1, 1, 1]
    for molecule_index, molecule in enumerate(group.sub_systems):
        np.testing.assert_array_equal(
            molecule.particle_indices, np.arange(molecule_index, 16, 2)
        )
        assert [subgroup.name for subgroup in molecule.sub_systems] == [
            f'group_{name}' for name in expected_names
        ]
        assert [len(subgroup.sub_systems) for subgroup in molecule.sub_systems] == (
            expected_counts
        )
        for subgroup, name in zip(molecule.sub_systems, expected_names, strict=True):
            assert {monomer.name for monomer in subgroup.sub_systems} == {name}
            assert set(subgroup.particle_indices) <= set(molecule.particle_indices)


def test_lammps_preserves_supplied_molecule_and_residue_type_names(monkeypatch):
    arrays, _ = interleaved_polymer_arrays()
    arrays['moltypes'][:] = 'PEO_chain'
    arrays['resnames'] = np.array(['CH3'] * 4 + ['EO'] * 8 + ['CH3'] * 4)
    writer = LammpsArchiveWriter()
    writer._mdanalysistraj_parser = SimpleNamespace(
        universe=SimpleNamespace(atoms=SimpleNamespace(moltypes=[], resnames=[])),
        get=lambda key, default=None: arrays if key == 'atoms_info' else default,
    )
    monkeypatch.setattr(
        writer, '_typed_bonds_by_particle', lambda universe, types: None
    )
    simulation = Simulation(model_system=[ModelSystem()])

    writer._parse_molecular_hierarchy(simulation)

    group = simulation.model_system[0].sub_systems[0]
    assert group.name == 'group_PEO_chain'
    for molecule in group.sub_systems:
        assert molecule.name == 'PEO_chain'
        assert [subgroup.name for subgroup in molecule.sub_systems] == [
            'group_CH3',
            'group_EO',
        ]
        assert [
            {monomer.name for monomer in subgroup.sub_systems}
            for subgroup in molecule.sub_systems
        ] == [{'CH3'}, {'EO'}]


def test_group_indices_by_value_groups_only_selected_indices_once():
    values = np.array(['B', 'A', 'B', 'A', 'B', 'A'])

    grouped = group_indices_by_value([5, 2, 3, 0], values)

    assert list(grouped) == ['A', 'B']
    assert grouped['A'].tolist() == [3, 5]
    assert grouped['B'].tolist() == [0, 2]


def test_conflicting_names_are_disambiguated_and_retained():
    repeated_large = system_node('PEO', 'molecule', range(8))
    same_signature = system_node('PEO', 'molecule', range(8, 16))
    repeated_small = system_node('PEO', 'monomer', range(3))
    root = system_node(
        'group_PEO',
        'molecule_group',
        range(16),
        [repeated_large, same_signature, repeated_small],
    )

    disambiguate_system_names([root])

    assert root['name'] == 'group_PEO'
    assert repeated_small['name'] == 'PEO_1'
    assert repeated_large['name'] == 'PEO_0'
    assert same_signature['name'] == repeated_large['name']
    assert root['sub_systems'] == [
        repeated_large,
        same_signature,
        repeated_small,
    ]
    assert repeated_small['particle_indices'].tolist() == [0, 1, 2]
    assert repeated_large['particle_indices'].tolist() == list(range(8))


def test_same_type_and_size_names_are_left_unchanged():
    first = system_node('PEO', 'molecule', range(3))
    second = system_node('PEO', 'molecule', range(3, 6))

    disambiguate_system_names([first, second])

    assert first['name'] == second['name'] == 'PEO'


def test_same_sized_groups_under_different_chains_get_distinct_names():
    chains = []
    for chain_index, chain_name in enumerate(['Protein_chain_A', 'Protein_chain_B']):
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

    disambiguate_system_names(chains)

    chain_a_group = chains[0]['sub_systems'][0]['sub_systems'][0]
    chain_b_group = chains[1]['sub_systems'][0]['sub_systems'][0]
    assert chain_a_group['name'] == 'group_SA_0'
    assert chain_b_group['name'] == 'group_SA_1'
    assert [group['name'] for group in chain_a_group['sub_systems']] == [
        'SA_0',
        'SA_0',
    ]
    assert [group['name'] for group in chain_b_group['sub_systems']] == [
        'SA_1',
        'SA_1',
    ]
    assert chain_a_group['particle_indices'].tolist() == list(range(6))
    assert chain_b_group['particle_indices'].tolist() == list(range(6, 12))


def test_disambiguated_names_avoid_existing_names():
    first = system_node('PEO', 'molecule', range(2))
    other_name = system_node('PEO_1', 'molecule', range(2, 4))
    conflicting = system_node('PEO', 'monomer', [0])

    disambiguate_system_names([first, other_name, conflicting])

    assert first['name'] == 'PEO_0'
    assert other_name['name'] == 'PEO_1'
    assert conflicting['name'] == 'PEO_2'
