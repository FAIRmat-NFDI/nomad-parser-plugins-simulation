import os
from typing import Any
import sys
import json

from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from nomad.datamodel.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import (
    MetainfoParser,
    TextParser,
    XMLParser,
)
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow.general import (
    DensityConvergenceTarget,
    EnergyConvergenceTarget,
    ForceConvergenceTarget,
    PotentialConvergenceTarget,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)
from structlog.stdlib import (
    BoundLogger,
)
from nomad_simulations.schema_packages.basis_set import (
    APWPlaneWaveBasisSet,
    BasisSetContainer,
    MuffinTinRegion,
)
from nomad_simulations.schema_packages.basis_set import (
    APWOrbital,
    APWLocalOrbital,
)
from nomad_simulations.schema_packages import basis_set
from nomad_simulation_parsers.parsers.utils.general import (
    calculate_band_gap_from_occupations,
    search_files,
)
from nomad_simulation_parsers.schema_packages import exciting

from .eigval_parser import EigvalFileParser
from .info_parser import InfoFileParser

convergence_threshold_mapping = {
    'x_exciting_effective_potential_convergence': {
        'class': PotentialConvergenceTarget,
        'threshold_type': 'rms',
    },
    'x_exciting_energy_convergence': {
        'class': EnergyConvergenceTarget,
        'threshold_type': 'absolute',
    },
    'x_exciting_charge_convergence': {
        'class': DensityConvergenceTarget,
        'threshold_type': 'absolute',
        'type': 'charge_abs',
    },
    'x_exciting_IBS_force_convergence': {
        'class': ForceConvergenceTarget,
        'threshold_type': 'absolute',
    },
}

XC_FUNCTIONAL_MAP = {
    2: {
        'key': 'LDA_PZ',
        'components': ['LDA_C_PZ', 'LDA_X'],
    },
    3: {
        'key': 'LDA_PW',
        'components': ['LDA_C_PW', 'LDA_X'],
    },
    4: {
        'key': 'XALPHA',
        'components': ['LDA_C_XALPHA'],
    },
    5: {
        'key': 'VBH',
        'components': ['LDA_C_VBH'],
    },
    20: {
        'key': 'PBE',
        'components': ['GGA_C_PBE', 'GGA_X_PBE'],
    },
    21: {
        'key': 'PBE_R',
        'components': ['GGA_C_PBE', 'GGA_X_PBE_R'],
    },
    22: {
        'key': 'PBESOL',
        'components': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
    },
    26: {
        'key': 'WC',
        'components': ['GGA_C_PBE', 'GGA_X_WC'],
    },
    30: {
        'key': 'AM05',
        'components': ['GGA_C_AM05', 'GGA_X_AM05'],
    },
    300: {
        'key': 'AC_PBE',
        'components': ['GGA_C_BGCP', 'GGA_X_PBE'],
    },
    406: {
        'key': 'PBE0',
        'components': ['HYB_GGA_XC_PBEH'],
    },
    408: {
        'key': 'HSE06',
        'components': ['HYB_GGA_XC_HSE06'],
    },
}


class ExcitingMetainfoParser(MetainfoParser):
    pass


class InfoParser(TextParser):
    def get_xc_functional_key(self, xc_type: int) -> str | None:
        functional = XC_FUNCTIONAL_MAP.get(xc_type)
        return functional['key'] if functional else None

    def get_xc_functional_components(self, xc_type: int) -> list[dict[str, Any]]:
        functional = XC_FUNCTIONAL_MAP.get(xc_type)
        if functional is None:
            return []
        return [{'libxc': name} for name in functional['components']]
       
    def get_forces(self, source: dict[str, Any]) -> dict[str, Any]:
        return dict(
            forces=source.get('forces'),
            n_points=len(source.get('forces', [])),
            rank=[3],
        )

    def get_configurations(self, root: dict[str, Any]) -> list[dict[str, Any]]:
        configurations = [
            root[key] for key in ['groundstate', 'hybrids'] if root.get(key)
        ]
        optimization = root.get('structure_optimization')
        if optimization:
            configurations.extend(optimization.get('optimization_step', []))
            configurations.append(optimization)
        # set final to last scf iteration
        for config in configurations:
            if not config.get('final') and config.get('scf_iteration'):
                config['final'] = config.get('scf_iteration')[-1]

        mapped_configurations = [
            # TODO do we really want to assign the initial atomic positions?
            self.get_atoms(config.get('atomic_positions', {}))
            for config in configurations
        ]
        merged = []
        for n, config in enumerate(mapped_configurations):
            entry = {**config, **configurations[n]}
            # Particle identity is frame-independent, so keep `atoms`
            # (-> particle_states) on the first (topology) frame only; later
            # optimization steps carry positions only. This avoids duplicating
            # per-atom identity per frame, which bloats the archive and the
            # Elasticsearch index doc (FAIRmat-NFDI/nomad-simulations#474).
            if n != 0:
                entry.pop('atoms', None)
            merged.append(entry)
        return merged

    def get_atoms(self, source: dict[str, Any]) -> dict[str, Any]:
        positions = source.get('positions')
        initial = self.data.get('initialization', {})
        lattice_vectors = initial.get('lattice_vectors')
        if positions is not None and source.get('positions_format') == 'lattice':
            positions = (
                np.dot(positions, lattice_vectors.magnitude) * lattice_vectors.units
            )
        if positions is None:
            positions = []
            for species in initial.get('species', []):
                positions_specie = species.get('positions')
                if species.get('positions_format') == 'lattice':
                    positions_specie = np.dot(
                        positions_specie, lattice_vectors.magnitude
                    )
                positions.extend(positions_specie)
            positions = positions * lattice_vectors.units

        atoms = []
        exclude = ['positions', 'positions_format', 'radial_points']
        for species in initial.get('species', []):
            atom = {k: v for k, v in species.items() if k not in exclude}
            atoms.extend([atom] * len(species.get('positions', [])))
        if not atoms:
             atoms = [dict(symbol=s) for s in source.get('symbols')]

        return dict(
            positions=(
                np.array(positions, dtype=np.float64) * ureg.bohr
                if not hasattr(positions, 'magnitude')
                else positions
            ),
            atoms=atoms,
            lattice_vectors=lattice_vectors,
            periodic_boundary_conditions=[True, True, True]
            if lattice_vectors is not None
            else None,
        )

    def get_geometry_convergence(
        self, source: dict[str, Any]
    ) -> list[ForceConvergenceTarget]:
        structure_optimization = source.get('structure_optimization')
        if structure_optimization is None:
            return []
        threshold = structure_optimization.get('force_target')
        if threshold is None:
            return []
        threshold = threshold.to('newton')
        convergence = [
            ForceConvergenceTarget(
                threshold=threshold,
                threshold_type='maximum',
            )
        ]
        # convergence.extend(self.get_single_point_convergence(source))
        return convergence

    def get_single_point_convergence(self, source: dict[str, Any]) -> list:
        groundstate = source.get('groundstate')
        if groundstate is None:
            return []
        scf_iterations = groundstate.get('scf_iteration', [])
        if not scf_iterations:
            return []
        last_iteration = scf_iterations[-1]

        convergence_targets = []
        for key_, info_ in convergence_threshold_mapping.items():
            quantity = last_iteration.get(key_, None)
            if quantity is None:
                continue
            # Keep pint quantity; metainfo handles conversion based on
            # target Quantity units.
            threshold_value = quantity[1]

            target_kwargs = {
                'threshold': threshold_value,
                'threshold_type': info_['threshold_type'],
            }
            if 'type' in info_:
                target_kwargs['type'] = info_['type']
            target = info_['class'](**target_kwargs)
            convergence_targets.append(target)
        return convergence_targets

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        scf_steps = source.get('groundstate', {}).get('scf_iteration', [])
        energies = []
        wall_times = []
        delta_energies = []
        delta_potential = []
        delta_charge = []
        delta_force = []

        def safe_append(source, value_name, unit_conversion, out):
            # Append values only if they exist for the current step
            value = source.get(value_name)
            if value is None:
                return
            out.append(value.to(unit_conversion)[0])

        for idx, step in enumerate(scf_steps):
            energies.append(step.get('energy_total').to('joule'))
            wall_times.append(step.get('time_physical').to('seconds').magnitude)
            safe_append(step, 'x_exciting_energy_convergence', 'joule', delta_energies)
            safe_append(
                step,
                'x_exciting_effective_potential_convergence',
                'joule',
                delta_potential,
            )
            safe_append(step, 'x_exciting_charge_convergence', 'coulomb', delta_charge)
            safe_append(step, 'x_exciting_IBS_force_convergence', 'newton', delta_force)
        durations = []
        # compute duration by subtracting previous step from cumulative time
        for idx, time in enumerate(wall_times):
            if idx == 0:
                duration = time
            else:
                duration = time - wall_times[idx - 1]
            durations.append(duration)
        out = {'energies_total': energies, 'durations': durations}
        for name, values in zip(
            [
                'delta_energies_total',
                'delta_potential_rms',
                'delta_charge_abs',
                'delta_force_abs',
            ],
            [delta_energies, delta_potential, delta_charge, delta_force],
        ):
            if len(values) > 0:
                out[name] = values
        return out

    def get_fermi_energy(self, source: dict[str, Any]) -> Any:
        """Resolve Fermi energy from INFO payload (legacy-equivalent source)."""
        if source is None:
            return None
        groundstate = source.get('groundstate') or {}

        # Prefer finalized value when available.
        final = groundstate.get('final') or {}
        fermi = final.get('x_exciting_fermi_energy')
        if fermi is not None:
            return fermi

        # Fallback to last SCF iteration.
        scf_iterations = groundstate.get('scf_iteration') or []
        if scf_iterations:
            fermi = scf_iterations[-1].get('x_exciting_fermi_energy')
            if fermi is not None:
                return fermi
        return None

    def get_is_spin_polarized(self, spin_treatment: str | None) -> bool:
        return {
            'spin-unpolarised': False,
            'spin-polarised': True,
        }.get(spin_treatment)

    def get_smearing_kind(self, smearing_kind: str | None) -> str | None:
        if smearing_kind is None:
            return None
        if isinstance(smearing_kind, (list, tuple)):
            smearing_kind = ' '.join(str(item) for item in smearing_kind)
        smearing_kind = ' '.join(smearing_kind.split())
        return {
            'Fermi-Dirac': 'Fermi-Dirac',
            'Gaussian': 'Gaussian',
            'Methfessel-Paxton order 1, Phys. Rev. B 40, 3616 (1989)': (
                'Methfessel-Paxton 1',
            ), 
            'Methfessel-Paxton order 2, Phys. Rev. B 40, 3616 (1989)': (
                'Methfessel-Paxton 2',
            ),
            'Extended linear tetrahedron method': 'Tetrahedra-extended',
            'Square-wave impulse': 'Square-wave impulse',
        }.get(smearing_kind)

    def get_k_mesh(self) -> list[dict[str, Any]]:
        initialization = self.data.get('initialization', {})
        grid = initialization.get('kpoint_grid')
        n_points = initialization.get('x_exciting_number_kpoints')
        return [
            {
                'grid': grid if grid is not None else None,
                'n_total_points': np.prod(grid) if grid is not None else None,
                'n_reduced_points': n_points if n_points is not None else None,
                'dimensionality': len(grid) if grid is not None else None,
            }
        ]
    def get_muffin_tin_species_scope(
        self,
        species_symbol: str | None,
    ) -> list[str]:
        if species_symbol is None:
            return []

        initialization = self.data.get('initialization', {})
        species_list = initialization.get('species', [])

        references = []
        atom_index = 0

        for species in species_list:
            symbol = species.get('symbol')
            positions = species.get('positions')
            n_atoms = len(positions) if positions is not None else 0

            if symbol == species_symbol:
                references.extend(
                    f'#/data/model_system/0/particle_states/{i}'
                    for i in range(atom_index, atom_index + n_atoms)
                )

            atom_index += n_atoms

        return references

class InputXMLParser(XMLParser):
    def get_xc_functionals(self, xc_funcs: dict[str, str]) -> list[dict[str, str]]:
        return [dict(libxc=val, type=key) for key, val in xc_funcs.items()]


class BandstructureXMLParser(XMLParser):
    n_spin = 1

    def get_k_path(self, source: dict[str, Any]) -> dict[str, Any]:
        def as_list(value: Any) -> list[Any]:
            if value is None:
                return []
            if isinstance(value, list):
                return value
            return [value]

        def parse_coord(coord: Any) -> np.ndarray | None:
            if coord is None:
                return None
            if isinstance(coord, str):
                coord = coord.split()
            try:
                parsed = np.array(coord, dtype=float)
            except Exception:
                return None
            if parsed.shape != (3,):
                return None
            return parsed

        bandstructure = source.get('bandstructure', {})

        # Prefer explicit point coordinates from the first band.
        points: list[np.ndarray] = []
        bands = as_list(bandstructure.get('band'))
        if bands:
            band_points = as_list(bands[0].get('point'))
            for point in band_points:
                coord = parse_coord(point.get('@coord'))
                if coord is not None:
                    points.append(coord)

        k_path: dict[str, Any] = {}
        if points:
            k_path['n_line_points'] = len(points)
            k_path['points'] = np.array(points, dtype=float)

        # Capture high-symmetry labels and vertices when available.
        high_symmetry_names: list[str] = []
        high_symmetry_values: list[np.ndarray] = []
        for vertex in as_list(bandstructure.get('vertex')):
            label = vertex.get('@label')
            coord = parse_coord(vertex.get('@coord'))
            if label is None or coord is None:
                continue
            high_symmetry_names.append(label)
            high_symmetry_values.append(coord)

        if high_symmetry_names and high_symmetry_values:
            k_path['high_symmetry_path_names'] = high_symmetry_names
            k_path['high_symmetry_path_values'] = np.array(high_symmetry_values)

        return k_path

    def get_bandstructures(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        # TODO determine format for spin pol case
        energies = [
            p['@eval'] for b in source['bandstructure']['band'] for p in b['point']
        ]
        n_spin = source.get('n_spin', self.n_spin)
        n_band = len(source['bandstructure']['band']) // n_spin
        n_kpoints = len(source['bandstructure']['band'][0]['point'])
        energies = np.array(energies, dtype=float).reshape((n_spin, n_band, n_kpoints))
        k_path = self.get_k_path(source)
        return [
            dict(
                energies=e.T * ureg.hartree,
                n_states=n_band,
                n_kpoints=n_kpoints,
                k_path=k_path,
            )
            for e in energies
        ]

    def reshape_coords(self, source: list[str]) -> np.ndarray:
        return np.array([v.split() for v in source], dtype=float)


class DosXMLParser(XMLParser):
    def to_float(self, source: list[str]) -> np.ndarray:
        return np.array(source, dtype=float)

    def get_dos(self, source: list[dict[str, Any]]) -> dict[str, Any]:
        return dict(
            dos=np.array([p['@dos'] for p in source.get('point', [])], dtype=float),
            energy=np.array([p['@e'] for p in source.get('point', [])], dtype=float),
        )


class EigvalParser(TextParser):
    def get_eigenvalues(self, source: dict[str, Any]):
        eigs_occs = source.get('eigenvalues_occupancies')
        eigs = np.array([v.get('eigenvalues') for v in eigs_occs])
        occs = np.array([v.get('occupancies') for v in eigs_occs])

        return [
            dict(
                eigenvalues=eigs[:, spin, :],
                occupancies=occs[:, spin, :],
                # n_states printed on file is actual n of states * n spin channels
                n_states=len(eigs[0][spin]),
            )
            for spin in range(len(eigs[0]))
        ]

    def get_band_gaps(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        """Calculate band gaps across all k-points for each spin channel."""
        eigs_occs = source.get('eigenvalues_occupancies') or []
        if not eigs_occs:
            return []

        n_spin = len(eigs_occs[0].get('eigenvalues', []))
        if n_spin == 0:
            return []

        band_gaps = []
        for spin in range(n_spin):
            # Concatenate over k-points so the gap is computed globally.
            eigenvalues = []
            occupations = []
            for kpoint in eigs_occs:
                eigs = np.asarray(kpoint.get('eigenvalues'))
                occs = np.asarray(kpoint.get('occupancies'))
                if eigs.size == 0 or occs.size == 0:
                    continue
                if spin >= eigs.shape[0] or spin >= occs.shape[0]:
                    continue
                eigenvalues.append(eigs[spin])
                occupations.append(occs[spin])
            if not eigenvalues:
                continue

            gap = calculate_band_gap_from_occupations(
                np.concatenate(eigenvalues), np.concatenate(occupations)
            )
            if gap is None:
                continue
            if n_spin > 1:
                gap['spin_channel'] = spin
            band_gaps.append(gap)

        return band_gaps

class SpeciesXMLParser(XMLParser):

    # def get_l_channels(self, basis: dict | None) -> list[dict]:
    #     if not basis:
    #         return []

    #     print('DEBUG basis:', basis)
    #     def as_list(value):
    #         if value is None:
    #             return []
    #         return value if isinstance(value, list) else [value]

    #     channels = {}

    #     for custom in as_list(basis.get('custom')):
    #         l = int(custom['@l'])

    #         channel = channels.setdefault(
    #             l, {'name': l, 'orbitals': []}
    #         )

    #         channel['orbitals'].append({
    #             'm_def': basis_set.APWOrbital.m_def.qualified_name(),
    #             **custom,
    #         })

    #     for lo in as_list(basis.get('lo')):
    #         l = int(lo['@l'])

    #         channel = channels.setdefault(
    #             l, {'name': l, 'orbitals': []}
    #         )

    #         channel['orbitals'].append({
    #             'm_def': basis_set.APWLocalOrbital.m_def.qualified_name(),
    #             **lo,
    #         })

    #     return [channels[l] for l in sorted(channels)]

    # def get_l_max(self, basis: dict | None) -> int | None:
    #     if not basis:
    #         return None

    #     l_values = []

    #     for key in ('custom', 'lo'):
    #         items = basis.get(key, [])

    #         if isinstance(items, dict):
    #             items = [items]

    #         for item in items:
    #             if item.get('@l') is not None:
    #                 l_values.append(int(item['@l']))

    #     return max(l_values) if l_values else None 
    def get_l_channels(self, basis: dict | None) -> list[dict]:
        if not basis:
            return []

        l_max = getattr(self, 'lmaxapw', None)
        if l_max is None:
            return []
        def as_list(value):
                if value is None:
                    return []
                return value if isinstance(value, list) else [value]

        default = basis.get('default')

        customs = {
            int(item['@l']): item
            for item in as_list(basis.get('custom'))
        }

        local_orbitals = {}
        for lo in as_list(basis.get('lo')):
            l = int(lo['@l'])
            local_orbitals.setdefault(l, []).append(lo)

        channels = []

        for l in range(int(l_max) + 1):
            orbitals = []

            # Use custom if present, otherwise default
            apw = customs.get(l, default)

            if apw is not None:
                orbitals.append({
                    '@type': 'lapw',
                    **apw,
                    '@l': str(l),
                })

            # Add all explicit local orbitals for this l
            orbitals.extend(local_orbitals.get(l, []))

            channels.append({
                'name': l,
                'orbitals': orbitals,
            })

        return channels
        # def as_list(value):
        #     if value is None:
        #         return []
        #     return value if isinstance(value, list) else [value]

        # channels = {}

        # for key in ('custom', 'lo'):
        #     for orbital in as_list(basis.get(key)):
        #         l = int(orbital['@l'])

        #         channel = channels.setdefault(
        #             l, {'name': l, 'orbitals': []}
        #         )

        #         channel['orbitals'].append(orbital)

        # return [channels[l] for l in sorted(channels)]
        
    def get_l_channel_name(self, l: str | None) -> int | None:
        return int(l) if l is not None else None

    # def get_apw_orbital_energy_parameter(
    #     self, trial_energy: str | None
    # ):
    #     if trial_energy is None:
    #         return None

    #     return np.array([float(trial_energy)]) * ureg.hartree

#     def get_apw_orbital_differential_order(
#     self, custom: dict | None
# ) -> np.ndarray | None:
#     if not custom:
#         return None

#     orbital_type = custom.get('@type')

#     orders = {
#         'apw': [0],
#         'lapw': [0, 1],
#         'slapw': [0, 1, 2],
#     }.get(orbital_type)

#     if orders is None:
#         return None

#     return np.array(orders, dtype=np.int32)


#     def get_apw_orbital_energy_parameter(
#         self, custom: dict | None
#     ):
#         if not custom or custom.get('@trialEnergy') is None:
#             return None

#         orders = self.get_apw_orbital_differential_order(custom)

#         if orders is None:
#             return None

#         return (
#             np.full(len(orders), float(custom['@trialEnergy']))
#             * ureg.hartree
#         )

    def get_local_orbital_energy_parameter(
        self, wf: dict | list[dict] | None
    ):
        if wf is None:
            return None

        if isinstance(wf, dict):
            wf = [wf]

        return np.array(
            [float(item['@trialEnergy']) for item in wf]
        ) * ureg.hartree

    def get_local_orbital_differential_order(
        self, wf: dict | list[dict] | None
    ):
        if wf is None:
            return None

        if isinstance(wf, dict):
            wf = [wf]

        return np.array(
            [int(item['@matchingOrder']) for item in wf],
            dtype=np.int32,
        )

    def get_local_orbital_energy_status(
        self, wf: dict | list[dict] | None
    ):
        mapping = {
            'false': 'fixed',
            'true': 'post-optimization',
        }
        if wf is None:
            return None

        if isinstance(wf, dict):
            wf = [wf]

        return [
                mapping.get(str(item.get('@searchE')))
                for item in wf
            ]

    def get_energy_status(self, orbital: dict | None):

        if not orbital:
            return None
        mapping = {
            'false': 'fixed',
            'true': 'post-optimization',
        }

        print('DEBUG energy_status 1:')
        if 'wf' in orbital:
            return self.get_local_orbital_energy_status(orbital['wf'])

        energy_status = orbital.get('@searchE')
        n_terms = {'apw': 1,'lapw': 2,}.get(orbital.get('@type'))
        status = mapping.get(str(orbital.get('@searchE')).lower())

        if status is None or n_terms is None:
            return None

        print('DEBUG energy_status:', energy_status)
        return [status] * n_terms

    def get_orbital_energy_parameter(self, orbital: dict | None):
        if not orbital:
            return None

        if 'wf' in orbital:
            return self.get_local_orbital_energy_parameter(orbital['wf'])

        energy = orbital.get('@trialEnergy')
        orbital_type = orbital.get('@type')

        n_terms = {'apw': 1, 'lapw': 2}.get(orbital_type)

        if energy is None or n_terms is None:
            return None

        return np.full(n_terms, float(energy)) * ureg.hartree


    def get_orbital_differential_order(self, orbital: dict | None):
        if not orbital:
            return None

        if 'wf' in orbital:
            return self.get_local_orbital_differential_order(orbital['wf'])

        orders = {
            'apw': [0],
            'lapw': [0, 1],
        }.get(orbital.get('@type'))

        return np.array(orders, dtype=np.int32) if orders is not None else None

    # def get_l_channels(self, basis):
    #     channels = {}

    #     for apw in basis.get('custom', []):
    #         l = int(apw['@l'])
    #         channels.setdefault(l, {'name': l, 'orbitals': []})

    #         channels[l]['orbitals'].append({
    #             'kind': 'apw',
    #             'type': apw['@type'],
    #             'trial_energy': float(apw['@trialEnergy']),
    #         })

    #     for lo in basis.get('lo', []):
    #         l = int(lo['@l'])
    #         channels.setdefault(l, {'name': l, 'orbitals': []})

    #         wavefunctions = lo['wf']
    #         if isinstance(wavefunctions, dict):
    #             wavefunctions = [wavefunctions]

   

class ExcitingArchiveWriter(ArchiveWriter):
    def write_to_archive(self) -> None:  # noqa: PLR0912, PLR0915
        maindir = os.path.dirname(self.mainfile)
        mainbase = os.path.basename(self.mainfile)

        # mainfile INFO.OUT parser
        info_parser = InfoParser(logger=self.logger, text_parser=InfoFileParser())
        info_parser.filepath = self.mainfile

        data_parser = ExcitingMetainfoParser(
            logger=self.logger, data_object=Simulation(program=Program(name='exciting'))
        )
        data_parser.annotation_key = exciting.INFO_KEY
        info_parser.convert(data_parser)
        species_list = info_parser.data.get('initialization', {}).get('species', [])

        for i, species in enumerate(species_list):
            print(f'DEBUG species {i}:', species.get('symbol'))
            print(
                f'DEBUG muffin-tin radius {i}:',
                species.get('muffin_tin_radius'),
            )   
        if not info_parser.data:
            info_out = os.path.join(maindir, 'INFO.OUT')

            if (
                os.path.isfile(info_out)
                and os.path.abspath(info_out)
                != os.path.abspath(self.mainfile)
            ):
                info_parser.filepath = info_out

                info_parser.convert(
                    data_parser,
                    update_mode='merge',
                )


        simulation = data_parser.data_object

        regions = [
            component
            for method in simulation.model_method or []
            for setting in method.numerical_settings or []
            if isinstance(setting, BasisSetContainer)
            for component in setting.basis_set_components or []
            if isinstance(component, MuffinTinRegion)
        ]

        species_list = info_parser.data.get('initialization', {}).get('species', [])

        if len(regions) != len(species_list):
            self.logger.warning('Muffin-tin regions and species count do not match.')
        else:
            for species, region in zip(species_list, regions):
                filename = species.get('file')

                if not filename:
                    continue

                xml_path = Path(maindir) / filename.strip()

                if not xml_path.is_file():
                    continue

                xml_parser = SpeciesXMLParser(filepath=str(xml_path))
                xml_parser.lmaxapw = info_parser.data[
                    'initialization'
                ]['x_exciting_lmaxapw']
                xml_target = ExcitingMetainfoParser(
                    logger=self.logger,
                    data_object=MuffinTinRegion(),
                )
                xml_target.annotation_key = exciting.SPECIES_XML_KEY
                xml_parser.convert(xml_target, update_mode='merge')
                xml_region = xml_target.data_object

                print(
                    'DEBUG XML l_channels:',
                    [(channel.name, channel.m_def.name) for channel in xml_region.l_channels],
                )
                try:
                    #xml_parser.convert(xml_target)
                    xml_region = xml_target.data_object

                    # l_max = xml_target.data_object.l_max

                    # if l_max is not None:
                    #     region.l_max = l_max
                    # xml_region = xml_target.data_object

                    if xml_region.l_max is not None:
                        region.l_max = xml_region.l_max

                    region.l_channels.clear()

                    for channel in xml_region.l_channels:
                        region.l_channels.append(channel)

                finally:
                    xml_parser.close()
                    xml_target.close()
        input_xml_files = (
            search_files('input.xml', maindir, re_pattern=mainbase)
            if not self.archive.m_xpath('data.model_method[0].xc.components')
            else []
        )
        if input_xml_files:
            input_xml_parser = InputXMLParser(
                filepath=input_xml_files[0], logger=self.logger
            )
            data_parser.annotation_key = exciting.INPUT_XML_KEY
            input_xml_parser.convert(data_parser)
            input_xml_parser.close()

        # eigenvalues from eigval.out
        eigval_files = search_files('EIGVAL.OUT', maindir, re_pattern=mainbase)
        if eigval_files:
            eigval_parser = EigvalParser(
                filepath=eigval_files[0],
                logger=self.logger,
                text_parser=EigvalFileParser(),
            )
            data_parser.annotation_key = exciting.EIGVAL_KEY
            eigval_parser.convert(data_parser, update_mode='merge')
            eigval_parser.close()

        # bandstructure from bandstructure.xml
        bandstructure_files = search_files(
            'bandstructure.xml', maindir, re_pattern=mainbase
        )
        if bandstructure_files:
            bandstructure_parser = BandstructureXMLParser(
                filepath=bandstructure_files[0], logger=self.logger
            )
            # TODO set n_spin from info
            data_parser.annotation_key = exciting.BANDSTRUCTURE_XML_KEY
            bandstructure_parser.convert(data_parser, update_mode='merge')
            bandstructure_parser.close()

        # dos from dos.xml
        dos_files = search_files('dos.xml', maindir, re_pattern=mainbase)
        if dos_files:
            dos_parser = DosXMLParser(filepath=dos_files[0], logger=self.logger)
            data_parser.annotation_key = exciting.DOS_XML_KEY
            dos_parser.convert(data_parser, update_mode='merge')
            dos_parser.close()

        self.archive.data = data_parser.data_object

        # Apply legacy-equivalent Fermi reference shift to electronic spectra.
        fermi_energy = info_parser.get_fermi_energy(info_parser.data)
        if fermi_energy is not None:
            for output in self.archive.data.outputs or []:
                for band_structure in output.electronic_band_structures or []:
                    if band_structure.value is not None:
                        band_structure.value = band_structure.value + fermi_energy
                    if band_structure.highest_occupied is None:
                        band_structure.highest_occupied = fermi_energy

                for dos in output.electronic_dos or []:
                    if dos.energies is not None and dos.energies.points is not None:
                        dos.energies.points = dos.energies.points + fermi_energy
                    if dos.energies_origin is None:
                        dos.energies_origin = fermi_energy

        # workflow section
        # populate geometry optimization if present
        if info_parser.text_parser.has_geometry_optimization():
            workflow = GeometryOptimization()
            workflow.method = GeometryOptimizationMethod()

            # Populate method using GEO_OPT_KEY annotations
            data_parser.data_object = workflow.method
            data_parser.annotation_key = exciting.GEO_OPT_KEY
            info_parser.convert(data_parser)

            # TODO: Investigate if mapping annotations can handle object
            # instantiation.
            # Currently, convergence targets are populated manually because:
            # 1. The mapping annotation system expects dict data from parsed files
            # 2. Our get_geometry_convergence() returns fully-formed metainfo objects
            # 3. The mapper cannot directly assign these objects; it tries to map
            # their fields.
            # Options to explore:
            # - Modify mapper to detect and handle object instances
            # - Change parser methods to return dicts that mapper can transform
            # - Keep manual population (current approach - clearer and more explicit)
            #
            # Keep this after info_parser.convert(); see the
            # `add_mapping_annotation` docstring for why the ordering matters.
            source_data = info_parser.data
            if source_data:
                convergence_targets = info_parser.get_geometry_convergence(source_data)
                if convergence_targets:
                    workflow.method.convergence_targets = convergence_targets

                # Also get single point convergence targets
                sp_convergence = info_parser.get_single_point_convergence(source_data)
                if sp_convergence:
                    workflow.method.single_point_convergence_targets = sp_convergence

            # Set the workflow with populated method
            workflow.method = data_parser.data_object
            self.archive.workflow2 = workflow
        else:  # here should come more standard workflows - for now only single point
            workflow = SinglePoint()
            workflow.method = SinglePointMethod()

            # Populate method using INFO_KEY annotations
            data_parser.data_object = workflow.method
            data_parser.annotation_key = exciting.INFO_KEY
            info_parser.convert(data_parser)

            # TODO: Same as above - manually populate convergence targets
            source_data = info_parser.data
            if source_data:
                sp_convergence = info_parser.get_single_point_convergence(source_data)
                if sp_convergence:
                    workflow.method.convergence_targets = sp_convergence

            # Set the workflow with populated method
            workflow.method = data_parser.data_object
            self.archive.workflow2 = workflow

        # close parsers
        info_parser.close()
        data_parser.close()


class ExcitingParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = ExcitingArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = None,
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
