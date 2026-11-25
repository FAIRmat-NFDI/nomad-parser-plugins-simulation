"""Build molecule and monomer hierarchies from particle arrays."""

from collections.abc import Callable, Iterable
from typing import Any

import numpy as np
from nomad_simulations.schema_packages.utils import get_composition


def group_indices_by_value(
    indices: Iterable[int], values: np.ndarray
) -> dict[Any, np.ndarray]:
    """Group a subset of particle indices by the corresponding array value."""
    indices = np.sort(np.asarray(indices))
    grouped_indices: dict[Any, list[int]] = {}
    for index in indices:
        grouped_indices.setdefault(values[index], []).append(index)

    return {
        value: np.asarray(grouped_indices[value], dtype=indices.dtype)
        for value in sorted(grouped_indices)
    }


def disambiguate_system_names(systems: Iterable[dict[str, Any]]) -> None:
    """Disambiguate names across hierarchy roles, particle counts, and parents.

    Equivalent occurrences retain the same name. Numeric suffixes distinguish
    conflicting hierarchy contexts.
    """
    systems = list(systems)
    reserved_names = set()
    names_by_context: dict[str, dict[tuple, list[dict[str, Any]]]] = {}

    def collect_contexts(
        siblings: Iterable[dict[str, Any]], parent_path: tuple = ()
    ) -> None:
        for system in siblings:
            base_name = system['name']
            reserved_names.add(base_name)
            signature = (
                system['branch_label'],
                len(system['particle_indices']),
                parent_path,
            )
            contexts = names_by_context.setdefault(base_name, {})
            contexts.setdefault(signature, []).append(system)
            child_path = parent_path + (
                (system['branch_label'], base_name, len(system['particle_indices'])),
            )
            collect_contexts(system.get('sub_systems', []), child_path)

    collect_contexts(systems)

    assigned_names = set()
    for base_name, contexts in names_by_context.items():
        if len(contexts) == 1:
            unique_names = [base_name]
        else:
            unique_names = []
            counter = 0
            for _ in contexts:
                unique_name = f'{base_name}_{counter}'
                while unique_name in reserved_names or unique_name in assigned_names:
                    counter += 1
                    unique_name = f'{base_name}_{counter}'
                unique_names.append(unique_name)
                counter += 1

        for matching_systems, unique_name in zip(contexts.values(), unique_names):
            assigned_names.add(unique_name)
            for system in matching_systems:
                system['name'] = unique_name


def infer_residue_types(
    particle_arrays: dict[str, np.ndarray],
    bonds_by_particle: list[list[tuple[int, str]]] | None,
) -> np.ndarray:
    """Infer missing residue type names from ordered, typed local topology.

    The caller must supply bonds in the same global particle-index ordering as
    the arrays. Internal connectivity and typed boundary bonds distinguish end
    groups from repeated interior residues. Particle ordering is part of the
    signature. If typed bonds or particle types are unavailable, preserve the
    original residue labels.
    """
    names = np.asarray(particle_arrays['resnames'], dtype=object).copy()
    types = particle_arrays['types']
    if bonds_by_particle is None or len(types) != len(names):
        return names

    type_by_signature: dict[tuple, str] = {}
    molecule_groups = group_indices_by_value(
        np.arange(len(names)), particle_arrays['moltypes']
    )
    for group_indices in molecule_groups.values():
        molecules = group_indices_by_value(group_indices, particle_arrays['molnums'])
        for molecule_indices in molecules.values():
            residues = group_indices_by_value(
                molecule_indices, particle_arrays['resids']
            )
            for indices in residues.values():
                local_indices = {int(index): slot for slot, index in enumerate(indices)}
                internal_bonds = []
                external_bonds = []
                for index, slot in local_indices.items():
                    for neighbor, bond_type in bonds_by_particle[index]:
                        neighbor_slot = local_indices.get(neighbor)
                        if neighbor_slot is None:
                            external_bonds.append(
                                (slot, str(types[neighbor]), bond_type)
                            )
                        elif slot < neighbor_slot:
                            internal_bonds.append((slot, neighbor_slot, bond_type))
                signature = (
                    tuple(str(types[index]) for index in indices),
                    tuple(sorted(internal_bonds)),
                    tuple(sorted(external_bonds)),
                )
                names[indices] = type_by_signature.setdefault(
                    signature, str(len(type_by_signature) + 1)
                )

    return names


def build_molecular_hierarchy(
    particle_arrays: dict[str, np.ndarray],
    chemical_formula: Callable[[np.ndarray], str] | None = None,
) -> dict[Any, dict[str, Any]]:
    """Build molecule/type groups with one section per physical occurrence.

    ``moltypes``, ``molnums``, ``resids``, and ``resnames`` must be aligned
    per-particle arrays. ``molnums`` and ``resids`` determine occurrence membership;
    names identify component types. All indices remain global, including interleaved
    memberships and residue IDs reused in different molecules. The optional
    formula callback supplies format-specific chemistry for terminal sections.
    """
    n_particles = len(particle_arrays['moltypes'])
    if any(
        len(particle_arrays[key]) != n_particles
        for key in ('molnums', 'resids', 'resnames')
    ):
        raise ValueError('Molecular hierarchy arrays must have equal particle counts.')

    def node(name: Any, role: str, indices: np.ndarray, **kwargs) -> dict[str, Any]:
        return {
            'name': str(name),
            'branch_label': role,
            'particle_indices': indices,
            'is_molecule': role == 'molecule',
            **kwargs,
        }

    hierarchy = {}
    molecule_groups = group_indices_by_value(
        np.arange(n_particles), particle_arrays['moltypes']
    )
    for moltype, group_indices in molecule_groups.items():
        molecules = group_indices_by_value(group_indices, particle_arrays['molnums'])
        molecule_group = node(
            f'group_{moltype}',
            'molecule_group',
            group_indices,
            sub_systems=[],
        )
        hierarchy[moltype] = molecule_group
        for molecule_indices in molecules.values():
            molecule = node(moltype, 'molecule', molecule_indices)
            molecule_group['sub_systems'].append(molecule)
            residues = group_indices_by_value(
                molecule_indices, particle_arrays['resids']
            )
            if len(residues) <= 1:
                if chemical_formula is not None:
                    molecule['composition_formula'] = chemical_formula(molecule_indices)
                continue

            residues_by_type: dict[Any, list[np.ndarray]] = {}
            for indices in residues.values():
                restype = particle_arrays['resnames'][indices[0]]
                residues_by_type.setdefault(restype, []).append(indices)
            molecule['sub_systems'] = []
            for restype in sorted(residues_by_type):
                typed_residues = residues_by_type[restype]
                monomer_group = node(
                    f'group_{restype}',
                    'monomer_group',
                    np.sort(np.concatenate(typed_residues)),
                    sub_systems=[],
                )
                molecule['sub_systems'].append(monomer_group)
                for indices in typed_residues:
                    monomer = node(restype, 'monomer', indices)
                    if chemical_formula is not None:
                        monomer['composition_formula'] = chemical_formula(indices)
                    monomer_group['sub_systems'].append(monomer)

    disambiguate_system_names(hierarchy.values())

    def set_parent_composition(system: dict[str, Any]) -> None:
        children = system.get('sub_systems', [])
        if not children:
            return
        # Construct structural formulas only after all type names are final.
        # Leaf chemistry supplied by the reader's callback is left untouched.
        system['composition_formula'] = get_composition(
            children_names=[child['name'] for child in children]
        )
        for child in children:
            set_parent_composition(child)

    for molecule_group in hierarchy.values():
        set_parent_composition(molecule_group)
    return hierarchy
