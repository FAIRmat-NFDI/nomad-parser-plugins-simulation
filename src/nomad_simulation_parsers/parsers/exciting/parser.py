import os
from typing import Any

import numpy as np
from nomad.datamodel.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import (
    MetainfoParser,
    TextParser,
    XMLParser,
)
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.workflow.general import (
    ChargeConvergenceTarget,
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

from nomad_simulation_parsers.parsers.utils.general import (
    calculate_band_gap_from_occupations,
    search_files,
)
from nomad_simulation_parsers.schema_packages import exciting

from .eigval_parser import EigvalFileParser
from .info_parser import InfoFileParser

LOGGER = get_logger(__name__)

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
        'class': ChargeConvergenceTarget,
        'threshold_type': 'absolute',
    },
    'x_exciting_IBS_force_convergence': {
        'class': ForceConvergenceTarget,
        'threshold_type': 'absolute',
    },
}


class ExcitingMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class InfoParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_xc_functionals(self, xc_type: int) -> list[dict[str, Any]]:
        xc_functional_map = {
            2: ['LDA_C_PZ', 'LDA_X_PZ'],
            3: ['LDA_C_PW', 'LDA_X_PZ'],
            4: ['LDA_C_XALPHA'],
            5: ['LDA_C_VBH'],
            20: ['GGA_C_PBE', 'GGA_X_PBE'],
            21: ['GGA_C_PBE', 'GGA_X_PBE_R'],
            22: ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
            26: ['GGA_C_PBE', 'GGA_X_WC'],
            30: ['GGA_C_AM05', 'GGA_C_AM05'],
            300: ['GGA_C_BGCP', 'GGA_X_PBE'],
            406: ['HYB_GGA_XC_PBEH'],
            408: ['HYB_GGA_XC_HSE03'],
        }
        return [dict(libxc=name) for name in xc_functional_map.get(xc_type, [])]

    def get_forces(self, source: dict[str, Any]) -> dict[str, Any]:
        strucopt = source.get('structure_optimization')
        return dict(
            forces=strucopt.get('forces'),
            n_points=len(strucopt.get('forces', [])),
            rank=[3],
        )

    def get_energies(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        """
        Extract total energies from all configurations:
        groundstate, hybrid, and optimization steps.
        """
        energies = []

        # Get energies from groundstate and hybrid calculations
        for key in ['groundstate', 'hybrid']:
            config = source.get(key)
            if config:
                # Check if there's a final section with energy_total
                final = config.get('final', {})
                if final.get('energy_total'):
                    energies.append({'energy_total': final['energy_total']})
                # Otherwise check if energy_total is directly in the config
                elif config.get('energy_total'):
                    energies.append({'energy_total': config['energy_total']})

        # Get energies from geometry optimization steps
        optimization = source.get('structure_optimization')
        if optimization:
            opt_steps = optimization.get('optimization_step', [])
            for step in opt_steps:
                if step.get('energy_total'):
                    energies.append({'energy_total': step['energy_total']})
            # Add final optimization energy
            if optimization.get('energy_total'):
                energies.append({'energy_total': optimization['energy_total']})

        return energies

    def get_configurations(self, root: dict[str, Any]) -> list[dict[str, Any]]:
        configurations = [
            root[key] for key in ['groundstate', 'hybrid'] if root.get(key)
        ]
        optimization = root.get('structure_optimization')
        if optimization:
            configurations.extend(optimization.get('optimization_step', []))
            configurations.append(optimization)
        mapped_configurations = [
            self.get_atoms(config['atomic_positions'])
            for config in configurations
            if config.get('atomic_positions')
        ]
        if mapped_configurations:
            return mapped_configurations

        # Fallback for minimal outputs where no explicit atomic_positions blocks
        # are present in groundstate/hybrid sections.
        if self.data.get('initialization'):
            return [self.get_atoms({})]
        return []

    def get_atoms(self, source: dict[str, Any]) -> dict[str, Any]:
        positions = source.get('positions')
        initial = self.data.get('initialization', {})
        lattice_vectors = initial.get('lattice_vectors')
        if positions is not None and source.get('positions_format') == 'lattice':
            positions = np.dot(positions, lattice_vectors.magnitude)
        if positions is None:
            positions = []
            for species in initial.get('species', []):
                positions_specie = species.get('positions')
                if species.get('positions_format') == 'lattice':
                    positions_specie = np.dot(positions_specie, lattice_vectors)
                positions.extend(positions_specie)
        atoms = []
        exclude = ['positions', 'positions_format', 'radial_points']
        for species in initial.get('species', []):
            atom = {k: v for k, v in species.items() if k not in exclude}
            atoms.extend([atom] * len(species.get('positions', [])))
        if not atoms:
            atoms = [dict(symbol=s) for s in source.get('symbols')]
        return dict(
            positions=np.array(positions, dtype=float),
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

            target = info_['class'](
                threshold=threshold_value,
                threshold_type=info_['threshold_type'],
            )
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
                'delta_density_rms',
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


class InputXMLParser(XMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_xc_functionals(self, xc_funcs: dict[str, str]) -> list[dict[str, str]]:
        return [dict(libxc=val, type=key) for key, val in xc_funcs.items()]


class BandstructureXMLParser(XMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

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
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_float(self, source: list[str]) -> np.ndarray:
        return np.array(source, dtype=float)

    def get_dos(self, source: list[dict[str, Any]]) -> dict[str, Any]:
        return dict(
            dos=np.array([p['@dos'] for p in source.get('point', [])], dtype=float),
            energy=np.array([p['@e'] for p in source.get('point', [])], dtype=float),
        )


class EigvalParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

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


class ExcitingArchiveWriter(ArchiveWriter):
    def write_to_archive(self) -> None:  # noqa: PLR0912, PLR0915
        maindir = os.path.dirname(self.mainfile)
        mainbase = os.path.basename(self.mainfile)

        # mainfile INFO.OUT parser
        info_parser = InfoParser(text_parser=InfoFileParser())
        info_parser.filepath = self.mainfile

        data_parser = ExcitingMetainfoParser(data_object=Simulation())
        data_parser.annotation_key = exciting.INFO_KEY

        info_parser.convert(data_parser)

        # Legacy exciting behavior uses INFO.OUT as the canonical source for
        # reference energy metadata. If the selected mainfile is GW_INFO.OUT
        # and parsing yields no payload, fallback to sibling INFO.OUT.
        if not info_parser.data:
            info_out = os.path.join(maindir, 'INFO.OUT')
            if os.path.isfile(info_out) and os.path.abspath(
                info_out
            ) != os.path.abspath(self.mainfile):
                info_parser.filepath = info_out
                info_parser.convert(data_parser)

        # read xc functionals from input.xml
        input_xml_files = (
            search_files('input.xml', maindir, re_pattern=mainbase)
            if not self.archive.m_xpath('data.model_method[0].xc_functionals')
            else []
        )
        if input_xml_files:
            input_xml_parser = InputXMLParser(filepath=input_xml_files[0])
            data_parser.annotation_key = exciting.INPUT_XML_KEY
            input_xml_parser.convert(data_parser)
            input_xml_parser.close()

        # eigenvalues from eigval.out
        eigval_files = search_files('EIGVAL.OUT', maindir, re_pattern=mainbase)
        if eigval_files:
            eigval_parser = EigvalParser(
                filepath=eigval_files[0], text_parser=EigvalFileParser()
            )
            data_parser.annotation_key = exciting.EIGVAL_KEY
            eigval_parser.convert(data_parser, update_mode='merge@-1')
            eigval_parser.close()

        # bandstructure from bandstructure.xml
        bandstructure_files = search_files(
            'bandstructure.xml', maindir, re_pattern=mainbase
        )
        if bandstructure_files:
            bandstructure_parser = BandstructureXMLParser(
                filepath=bandstructure_files[0]
            )
            # TODO set n_spin from info
            data_parser.annotation_key = exciting.BANDSTRUCTURE_XML_KEY
            bandstructure_parser.convert(data_parser, update_mode='merge@-1')
            bandstructure_parser.close()

        # dos from dos.xml
        dos_files = search_files('dos.xml', maindir, re_pattern=mainbase)
        if dos_files:
            dos_parser = DosXMLParser(filepath=dos_files[0])
            data_parser.annotation_key = exciting.DOS_XML_KEY
            dos_parser.convert(data_parser, update_mode='merge@-1')
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
