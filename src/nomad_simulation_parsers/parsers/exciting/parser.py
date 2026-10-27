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
        # atoms = []
        # exclude = ['positions', 'positions_format', 'radial_points']

        # species_xml_by_number = getattr(self, 'species_xml_by_number', {})

        # for specie in initial.get('species', []):
        #     atom = {k: v for k, v in specie.items() if k not in exclude}

        #     xml_data = species_xml_by_number.get(int(specie['number']))
        #     if xml_data:
        #         atom.update({
        #             'atomic_number': xml_data['atomic_number'],
        #             'chemical_symbol': xml_data['chemical_symbol'],
        #             'mass': xml_data['mass'],
        #         })

        # atoms.extend([atom.copy() for _ in specie.get('positions', [])])

        atoms = []
        exclude = ['positions', 'positions_format', 'radial_points']
        for species in initial.get('species', []):
            atom = {k: v for k, v in species.items() if k not in exclude}
            atoms.extend([atom] * len(species.get('positions', [])))
        if not atoms:
             atoms = [dict(symbol=s) for s in source.get('symbols')]
        # print('DEBUG get_atoms called')
        # print('DEBUG XML species:', species_xml_by_number)
        # print('DEBUG INFO species:', initial.get('species', []))
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

class SpeciesXMLParser:
    def parse(self, filepath: Path) -> dict:
        root = ET.parse(filepath).getroot()
        sp = root.find('sp')

        if sp is None:
            return {}

        muffin_tin = sp.find('muffinTin')

        return {
            'atomic_number': abs(int(float(sp.get('z')))),
            'chemical_symbol': sp.get('chemicalSymbol'),
            'mass': float(sp.get('mass')) * ureg.electron_mass,
            'muffin_tin': {
                'radius': float(muffin_tin.get('radius')) * ureg.bohr,
            } if muffin_tin is not None else {},
        }

    def parse_all(
            self,
            calculation_dir: Path,
            species: list,
        ) -> dict[int, dict]:
            calculation_dir = Path(calculation_dir).resolve()
            species_data = {}

            for specie in species:
                number = int(specie['number'])
                symbol = specie['symbol']
                filename = specie['file'].strip()

                filepath = (calculation_dir / filename).resolve()

                if not filepath.is_relative_to(calculation_dir):
                    raise ValueError(f'Invalid species XML path: {filename}')

                if not filepath.is_file():
                    raise FileNotFoundError(
                        f'Species XML file not found: {filepath}'
                    )

                data = self.parse(filepath)

                if data['chemical_symbol'] != symbol:
                    raise ValueError(
                        f'Species mismatch in {filename}: '
                        f'expected {symbol}, '
                        f"found {data['chemical_symbol']}"
                    )

                if number in species_data:
                    raise ValueError(
                        f'Duplicate species number: {number}'
                    )

                species_data[number] = data

                print(
                    f'Species {number}: {symbol}, '
                    f'file={filepath}, '
                    f"MT radius={data['muffin_tin'].get('radius')}"
                )

            return species_data
        # for specie in sorted(species, key=lambda s: int(s['number'])):
        #     filename = specie['file'].strip()
        #     filepath = calculation_dir / filename

        #     if not filepath.is_file():
        #         species_data.append(None)
        #         continue

        #     data = self.parse(filepath)

        #     if data['chemical_symbol'] != specie['symbol']:
        #         raise ValueError(
        #             f'Species mismatch in {filename}: '
        #             f"expected {specie['symbol']}, "
        #             f"found {data['chemical_symbol']}"
        #         )

        #     species_data.append(data)
        # for specie in species:
        #     print(
        #         'Species:', specie['number'],
        #         'Symbol:', specie['symbol'],
        #         'File:', specie['file'],
        #     )
        # return species_data


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

        
        # ---------------------------------------------------------
        # 1. Convert INFO.OUT
        # ---------------------------------------------------------
        data_parser.annotation_key = exciting.INFO_KEY
        info_parser.convert(data_parser)
        
    

        simulation = data_parser.data_object

        print('\nDEBUG MODEL SYSTEM STRUCTURE START', file=sys.stderr)

        for i, model_system in enumerate(simulation.model_system or []):
            print(f'\nMODEL SYSTEM [{i}]', file=sys.stderr)

            print(
                json.dumps(
                    model_system.m_to_dict(),
                    indent=2,
                    default=str,
                ),
                file=sys.stderr,
            )

        print('\nDEBUG MODEL SYSTEM STRUCTURE END', file=sys.stderr)

        # Legacy fallback: if mainfile is GW_INFO.OUT and contains
        # no data, try the sibling INFO.OUT file.
        if not info_parser.data:
            info_out = os.path.join(maindir, 'INFO.OUT')

            if (
                os.path.isfile(info_out)
                and os.path.abspath(info_out)
                != os.path.abspath(self.mainfile)
            ):
                info_parser.filepath = info_out

                data_parser.annotation_key = exciting.INFO_KEY
                info_parser.convert(
                    data_parser,
                    update_mode='merge',
                )
         # ---------------------------------------------------------
        # 2. Read species XML files
        # ---------------------------------------------------------
        initial = info_parser.data.get('initialization', {})
        species = initial.get('species', [])
        print(
            f'DEBUG number of species: {len(species)}',
            file=sys.stderr,
        )
        for specie in species:
            symbol = specie.get('symbol', '?')
            filename = specie.get('file', '').strip()

            if not filename:
                print(
                    f'DEBUG missing filename for species {symbol}',
                    file=sys.stderr,
                )
                continue

            xml_path = Path(maindir) / filename

            if not xml_path.is_file():
                print(
                    f'DEBUG XML file not found: {xml_path}',
                    file=sys.stderr,
                )
                continue

            xml_parser = XMLParser()
            xml_parser.filepath = str(xml_path)

            # -----------------------------------------------------
            # 3. Inspect XML data before conversion
            # -----------------------------------------------------
            xml_data = xml_parser.data

            muffin_tin = (
                xml_data
                .get('spdb', {})
                .get('sp', {})
                .get('muffinTin', {})
            )

            radius_raw = muffin_tin.get('@radius')

            print(
                f'\nDEBUG BEFORE XML CONVERSION',
                file=sys.stderr,
            )
            print(
                f'  species: {symbol}',
                file=sys.stderr,
            )
            print(
                f'  file: {xml_path}',
                file=sys.stderr,
            )
            print(
                f'  radius raw: {radius_raw}',
                file=sys.stderr,
            )

            if radius_raw is None:
                print(
                    f'  WARNING: muffinTin radius missing for {symbol}',
                    file=sys.stderr,
                )
                continue

            # -----------------------------------------------------
            # 4. Convert XML using SPECIES_XML_KEY
            # -----------------------------------------------------
            data_parser.annotation_key = exciting.SPECIES_XML_KEY

            xml_parser.convert(
                data_parser,
                update_mode='merge',
            )
            radius = float(radius_raw) * ureg.bohr

            simulation = data_parser.data_object

            for method in simulation.model_method or []:
                for setting in method.numerical_settings or []:
                    if not isinstance(setting, BasisSetContainer):
                        continue

                    # Seleziona il contenitore che contiene APWPlaneWaveBasisSet
                    if not any(
                        isinstance(comp, APWPlaneWaveBasisSet)
                        for comp in setting.basis_set_components
                    ):
                        continue

                    # muffin_tin_region = MuffinTinRegion(
                    #     radius=radius,
                    # )

                    # setting.basis_set_components.append(muffin_tin_region)

                    # Recupera gli AtomsState già presenti nella simulazione
                    matching_atoms_states = []

                    for model_system in data_parser.data_object.model_system or []:
                        for atoms_state in model_system.particle_states or []:
                            if atoms_state.chemical_symbol == 'C':
                                matching_atoms_states.append(atoms_state)

                    # Crea il MuffinTinRegion con il raggio e i riferimenti agli atomi
                    muffin_tin_region = MuffinTinRegion(
                        radius=radius,
                        species_scope=matching_atoms_states,
                    )

                    # Aggiunge il componente al BasisSetContainer esistente
                    setting.basis_set_components.append(muffin_tin_region)

                    print(
                        f'DEBUG ADDED MuffinTinRegion: '
                        f'species={"C"}, '
                        f'radius={radius}, '
                        f'number of linked atoms={len(matching_atoms_states)}',
                        file=sys.stderr,
                    )
                    print(
                        f'DEBUG ADDED MuffinTinRegion: '
                        f'species={symbol}, radius={muffin_tin_region.radius}',
                        file=sys.stderr,
                    )

                    break
            print(
                f'DEBUG AFTER XML CONVERSION: {symbol}',
                file=sys.stderr,
            )

            # -----------------------------------------------------
            # 5. Inspect resulting BasisSetContainer objects
            # -----------------------------------------------------
            simulation = data_parser.data_object

            for method_index, method in enumerate(
                simulation.model_method or []
            ):
                numerical_settings = getattr(
                    method,
                    'numerical_settings',
                    [],
                ) or []

                for basis_index, basis in enumerate(
                    numerical_settings
                ):
                    if not hasattr(basis, 'basis_set_components'):
                        continue

                    print(
                        f'  BASIS [{method_index}, {basis_index}]: '
                        f'name={basis.name}',
                        file=sys.stderr,
                    )

                    components = basis.basis_set_components or []

                    print(
                        f'  number of components: {len(components)}',
                        file=sys.stderr,
                    )

                    for component in components:
                        if isinstance(
                            component,
                            APWPlaneWaveBasisSet,
                        ):
                            print(
                                '    APWPlaneWaveBasisSet: '
                                f'cutoff_fractional='
                                f'{component.cutoff_fractional}',
                                file=sys.stderr,
                            )

                        elif isinstance(
                            component,
                            MuffinTinRegion,
                        ):
                            print(
                                '    MuffinTinRegion: '
                                f'radius={component.radius}',
                                file=sys.stderr,
                            )

                        else:
                            print(
                                '    OTHER COMPONENT: '
                                f'{component.m_def.name}',
                                file=sys.stderr,
                            )

            print(
                'DEBUG END SPECIES\n',
                file=sys.stderr,
            )

        # for specie in species:
        #     xml_path = Path(maindir) / specie['file'].strip()

        #     if not xml_path.is_file():
        #         continue

        #     xml_parser = XMLParser()
        #     xml_parser.filepath = str(xml_path)

        # print('DEBUG species:', specie['symbol'])
        # print('DEBUG XML data:', xml_parser.data)
        # print(xml_parser.data['spdb']['sp']['muffinTin']['@radius'])
        # # # 2. Load species XML data
        # # species_xml_data = SpeciesXMLParser().parse_all(
        # #     calculation_dir=Path(maindir),
        # #     species=species,
        # # )

        # # # 4. Convert INFO.OUT after attaching XML data
        # data_parser.annotation_key = exciting.INFO_KEY
        # info_parser.convert(data_parser)

        # radius = xml_parser.data['spdb']['sp']['muffinTin']['@radius']

        # print('DEBUG radius:', radius, file=sys.stderr)

        # data_parser.annotation_key = exciting.SPECIES_XML_KEY

        # xml_parser.convert(
        #     data_parser,
        #     update_mode='merge',
        # )

        # # Parse species XML
        # species_xml_data = SpeciesXMLParser().parse_all(    
        #     calculation_dir=Path(maindir),
        #     species=species,)

        # species_xml_mapping_data = {
        #     'muffin_tin_regions': [
        #         {
        #             'species_number': number,
        #             'chemical_symbol': data['chemical_symbol'],
        #             'radius': data['muffin_tin']['radius'],
        #         }
        #         for number, data in sorted(species_xml_data.items())
        #         if data.get('muffin_tin')
        #     ]
        # }

        # print('DEBUG species XML mapping:', species_xml_mapping_data)
        # # Map species XML into MuffinTinRegion
        # data_parser.annotation_key = exciting.SPECIES_XML_KEY
        # print('DEBUG annotation key:', data_parser.annotation_key)
        # print('DEBUG XML data:', species_xml_mapping_data)
        # Legacy exciting behavior uses INFO.OUT as the canonical source for
        # reference energy metadata. If the selected mainfile is GW_INFO.OUT
        # and parsing yields no payload, fallback to sibling INFO.OUT.
        # if not info_parser.data:
        #     info_out = os.path.join(maindir, 'INFO.OUT')
        #     if os.path.isfile(info_out) and os.path.abspath(
        #         info_out
        #     ) != os.path.abspath(self.mainfile):
        #         info_parser.filepath = info_out
        #         info_parser.convert(data_parser, update_mode='merge')
        
        # Read species XML files
        # species = info_parser.data.get('initialization', {}).get('species', [])
        # species_xml_data = SpeciesXMLParser().parse_all(
        #     calculation_dir=maindir,
        #     species=species,
        # )
        # info_parser.species_xml_data = species_xml_data
        
        # read xc functionals from input.xml only if INFO.OUT did not already
        # populate them
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
