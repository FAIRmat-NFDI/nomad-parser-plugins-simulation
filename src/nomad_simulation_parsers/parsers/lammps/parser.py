import os
import re
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad_simulations.schema_packages.force_field import ForceField
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_method import ModelMethod
from nomad_simulations.schema_packages.model_system import ModelSystem
from nomad_simulations.schema_packages.outputs import TrajectoryOutputs, WorkflowOutputs
from nomad_simulations.schema_packages.properties.energies import (
    BaseEnergy,
    KineticEnergy,
    PotentialEnergy,
    TotalEnergy,
)
from nomad_simulations.schema_packages.properties.forces import TotalForce
from nomad_simulations.schema_packages.properties.thermodynamics import (
    Enthalpy,
    MassDensity,
    Pressure,
    Temperature,
    Volume,
)
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    MolecularDynamics,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.lammps.file_parsers import (
    DataParser,
    LogParser,
    get_unit,
)
from nomad_simulation_parsers.parsers.lammps.force_field import (
    attach_particle_parameters,
    build_force_field,
)
from nomad_simulation_parsers.parsers.lammps.md_workflow_method import (
    build_md_method,
    selected_run,
)
from nomad_simulation_parsers.parsers.lammps.trajectory_parsers import (
    TrajParser,
    TrajParsers,
    XYZTrajParser,
)
from nomad_simulation_parsers.parsers.utils.mdanalysisparser import MDAnalysisParser
from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser

_ENERGY_CONTRIBUTION_TYPES = {
    'e_pair': 'pair',
    'e_bond': 'bond',
    'e_angle': 'angle',
    'e_dihed': 'dihedral',
    'e_impro': 'improper',
    'e_coul': 'coulomb',
    'e_vdwl': 'van der Waals',
    'e_mol': 'molecular',
    'e_long': 'kspace long range',
    'e_tail': 'van der Waals long range',
}


class LammpsArchiveWriter(MDParser):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._log_parser = LogParser()
        self._aux_log_parser = LogParser()
        self._traj_parser = TrajParser()
        self._xyztraj_parser = XYZTrajParser()
        self._mdanalysistraj_parser = MDAnalysisParser(
            topology_format='DATA', format='LAMMPSDUMP'
        )
        self._data_parser = DataParser()
        self._resolved_unit_style: str | None = None
        self._bond_list: list[tuple[int, int]] | None = None
        self._trajectory_frame_count = 0

    def _require_physical_unit_style(self) -> str:
        """Resolve units from commands or same-log runtime evidence.

        A neighboring input is not sufficient evidence that it produced this
        log. Different command and runtime styles are also unsafe to combine.
        """
        unit_reader = LogParser()
        unit_reader.mainfile = self.mainfile
        try:
            declared_units = unit_reader.get('units') or []
        finally:
            unit_reader.close()

        if isinstance(declared_units, str):
            declared_units = [declared_units]
        runtime_units = []
        runtime_pattern = re.compile(r'^\s*Unit\s+style\s*:\s*([A-Za-z]+)\b', re.I)
        with open(self.mainfile, encoding='utf-8', errors='replace') as source:
            for line in source:
                match = runtime_pattern.match(line)
                if match:
                    runtime_units.append(match.group(1))

        if not declared_units and not runtime_units:
            raise ValueError(
                'LAMMPS unit style is unresolved in the log. Reduced LJ units '
                'are the LAMMPS default, and this parser cannot assign physical '
                'units without a supported command or same-log runtime unit style.'
            )

        unit_styles = {style.lower() for style in [*declared_units, *runtime_units]}
        if len(unit_styles) != 1:
            raise ValueError(
                f'Conflicting LAMMPS unit styles in the log: {sorted(unit_styles)}.'
            )

        unit_style = unit_styles.pop()
        if unit_style == 'lj':
            raise ValueError(
                'LAMMPS units lj are reduced and cannot be mapped to NOMAD '
                'physical units without system-specific scales. Parsing stopped '
                'before writing a simulation archive.'
            )

        physical_styles = {'real', 'metal', 'si', 'cgs', 'electron', 'micro', 'nano'}
        if unit_style not in physical_styles:
            raise ValueError(f'Unsupported LAMMPS unit style: {unit_style}.')
        return unit_style

    def apply_unit(self, value: Any, unit: str) -> float:
        if not hasattr(value, 'units'):
            value = value * self._log_parser.units.get(unit, 1)
        return value

    def parse_method(self, simulation: Simulation) -> None:
        """Attach source-backed force-field terms and prepare trajectory topology."""

        if self.traj_parsers[0].mainfile is None:
            return

        if self.traj_parsers.eval('n_frames') is None:
            return

        masses_data = self._data_parser.get('Masses', None)
        # Extract array from DataParser format: [(None, np.ndarray)]
        if masses_data and isinstance(masses_data, list) and len(masses_data) > 0:
            masses = (
                masses_data[0][1] if isinstance(masses_data[0], tuple) else masses_data
            )
        else:
            masses = None

        # Set masses on all trajectory parsers since eval() can return from any parser
        for parser in self.traj_parsers._parsers:
            if isinstance(parser, TrajParser):
                parser.masses = masses

        # TODO: find best place for first attempt to set _bond_list
        # Extract bond list from MDAnalysis universe if available
        if (
            self._mdanalysistraj_parser.mainfile is not None
            and self._mdanalysistraj_parser.universe is not None
        ):
            self._bond_list = [
                tuple(interaction['atom_indices'])
                for interaction in self._mdanalysistraj_parser.get_interactions()
                if interaction['type'] == 'bond'
            ]
        else:
            self._bond_list = None

        atoms_info = self.traj_parsers.eval('atoms_info')
        if isinstance(atoms_info, list):
            atoms_info = atoms_info[0] if atoms_info else None
        atom_ids = atoms_info.get('id') if isinstance(atoms_info, dict) else None
        method = build_force_field(
            self._log_parser,
            self._data_parser,
            atom_ids,
            get_unit(self._resolved_unit_style),
            self.logger,
        )
        simulation.model_method.append(method or ModelMethod())

    def _source_trajectory_step(self, frame_index: int) -> int | None:
        """Return a LAMMPS step ID when the frame reader exposes one."""
        if len(self.traj_parsers._parsers) != 1:
            return None
        reader = self.traj_parsers._parsers[0]
        if isinstance(reader, MDAnalysisParser) and (
            reader.options.get('format') != 'LAMMPSDUMP'
        ):
            return None
        raw_step = reader.get_step(frame_index)
        try:
            step = float(raw_step)
        except (TypeError, ValueError):
            return None
        if not np.isfinite(step) or not step.is_integer():
            return None
        if step < 0 or step > np.iinfo(np.int32).max:
            return None
        return int(step)

    def _validate_trajectory_data(self) -> bool:
        """Validate that trajectory data is available and extract basic info."""
        n_traj = self.traj_parsers.eval('n_frames')
        if n_traj is None:
            return False

        self._trajectory_frame_count = n_traj
        self.n_particles = [
            self.traj_parsers.eval('get_n_atoms', n) for n in range(n_traj)
        ]
        self.trajectory_steps = [
            step
            for n in range(n_traj)
            if (step := self.traj_parsers.eval('get_step', n)) is not None
        ]
        # MDParser deduplicates step IDs for its output lookup. A dump can
        # legitimately contain multiple frames at the same step after a reset;
        # sampling must use the physical frame count rather than unique steps.
        self.info['n_frames'] = n_traj
        return True

    def _parse_trajectory_frames(self, simulation: Simulation) -> None:
        """Parse all trajectory frames and create model systems."""

        def _get_quantity_with_units(
            traj_n: int, quantity: str, unit: str
        ) -> np.ndarray | None:
            """Get lattice vectors, velocities with units applied."""
            data = self.traj_parsers.eval(f'get_{quantity}', traj_n)
            if data is not None:
                return self.apply_unit(data, unit)
            return None

        def _extract_bond_list() -> None:
            """Extract bond list from data parser if not already set."""
            if self._bond_list is not None:
                return

            bonds = self._data_parser.get('Bonds', None)
            if bonds is None or bonds[0][1].size == 0:
                self._bond_list = None
                return

            # LAMMPS Bonds contain atom IDs, not zero-based indices into the
            # sorted trajectory arrays. Do not guess a mapping when IDs are
            # absent or inconsistent with the selected trajectory reader.
            atoms_info = self.traj_parsers.eval('atoms_info')
            if isinstance(atoms_info, list):
                atoms_info = atoms_info[0] if atoms_info else None
            raw_ids = atoms_info.get('id') if isinstance(atoms_info, dict) else None
            if raw_ids is None:
                self.logger.warning(
                    'Cannot map LAMMPS bonds without trajectory atom IDs.'
                )
                return

            atom_ids = np.asarray(raw_ids, dtype=float)
            bond_ids = np.asarray(bonds[0][1][:, 2:4], dtype=float)
            if (
                len(atom_ids) != self.traj_parsers.eval('get_n_atoms', 0)
                or not np.isfinite(atom_ids).all()
                or not np.isfinite(bond_ids).all()
                or not np.array_equal(atom_ids, atom_ids.astype(int))
                or not np.array_equal(bond_ids, bond_ids.astype(int))
            ):
                self.logger.warning('Cannot map inconsistent LAMMPS bond atom IDs.')
                return

            atom_id_to_index = {
                int(atom_id): index for index, atom_id in enumerate(atom_ids)
            }
            if len(atom_id_to_index) != len(atom_ids):
                self.logger.warning('Cannot map duplicate LAMMPS trajectory atom IDs.')
                return
            try:
                self._bond_list = [
                    (atom_id_to_index[int(a)], atom_id_to_index[int(b)])
                    for a, b in bond_ids
                ]
            except KeyError:
                self.logger.warning(
                    'Cannot map LAMMPS bond IDs missing from the trajectory.'
                )

        for traj_n in range(
            0, self._trajectory_frame_count, self.archive_sampling_rate
        ):
            # Extract and apply units to trajectory data
            lattice_vectors = _get_quantity_with_units(
                traj_n, 'lattice_vectors', 'distance'
            )
            velocities = _get_quantity_with_units(traj_n, 'velocities', 'velocity')
            # The 'dimension' command "must be used" to specify 2D simulations (https://docs.lammps.org/Howto_2d.html).
            # LAMMPS default is 3D.
            dimension = self._log_parser.get('dimension', 3)

            # Extract bond list for first frame only
            # TODO: add link to this in other frames
            if traj_n == 0:
                _extract_bond_list()

            particles_dict = {
                'lattice_vectors': lattice_vectors,
                'periodic_boundary_conditions': self.traj_parsers.eval(
                    'get_pbc', traj_n
                ),
                # Particle identity is identical in every frame; attach labels only to
                # the first (topology) frame so `particle_states` is stored and indexed
                # once, not per frame (FAIRmat-NFDI/nomad-simulations#474).
                'labels': self.traj_parsers.eval('get_atom_labels', traj_n)
                if traj_n == 0
                else None,
                'n_particles': self.traj_parsers.eval('get_n_atoms', traj_n),
                'positions': self.apply_unit(
                    self.traj_parsers.eval('get_positions', traj_n), 'distance'
                ),
                'velocities': velocities,
                'bond_list': self._bond_list if self._bond_list else None,
                'dimensions': dimension,
            }
            source_step = self._source_trajectory_step(traj_n)
            model_system = (
                ModelSystem(time_step=source_step) if source_step is not None else None
            )
            self.parse_trajectory_step(
                particles_dict, simulation, model_system=model_system
            )

    def _create_system_node(
        self, name: str | int, branch_label: str, particle_indices: np.ndarray, **kwargs
    ) -> ModelSystem:
        """
        Create a ModelSystem node with common setup.

        Args:
            name: System name
            branch_label: Hierarchy level label
            particle_indices: Indices of particles in this system
            **kwargs: Additional attributes: composition_formula, is_representative, ...
        """
        system = ModelSystem()
        system.name = str(name)
        system.branch_label = branch_label
        system.particle_indices = particle_indices

        # Set any additional attributes
        for key, value in kwargs.items():
            setattr(system, key, value)

        return system

    def _create_molecule(
        self, molecule: int, i_molecule: int, particle_arrays: dict
    ) -> ModelSystem:
        """Create a single molecule with its residues."""

        def _create_residue(
            res_id: int,
            restype: str,
            parent_system: ModelSystem,
            particle_arrays: dict,
        ) -> ModelSystem:
            """Create a single residue."""
            particle_indices = np.where(particle_arrays['resids'] == res_id)[0]
            particle_indices = np.intersect1d(
                particle_indices, parent_system.particle_indices
            )

            return self._create_system_node(
                name=restype,
                branch_label='monomer',
                particle_indices=particle_indices,
            )

        def _create_monomer_group(
            restype: str, parent_system: ModelSystem, particle_arrays: dict
        ) -> ModelSystem:
            """Create a monomer group with its constituent residues."""

            restype_indices = np.where(particle_arrays['resnames'] == restype)[0]
            particle_indices = np.intersect1d(
                restype_indices, parent_system.particle_indices
            )

            monomer_group = self._create_system_node(
                name=f'group_{restype}',
                branch_label='monomer_group',
                particle_indices=particle_indices,
            )

            # Add individual residues
            restype_resids = np.unique(
                particle_arrays['resids'][monomer_group.particle_indices]
            )
            for res_id in restype_resids:
                residue = _create_residue(
                    res_id, restype, monomer_group, particle_arrays
                )
                monomer_group.sub_systems.append(residue)

            return monomer_group

        def _add_residue_hierarchy(
            sec_molecule: ModelSystem, particle_arrays: dict
        ) -> None:
            """Add residue/monomer hierarchy to a molecule."""
            mol_resnames = particle_arrays['resnames'][sec_molecule.particle_indices]
            restypes = np.unique(mol_resnames)

            for restype in restypes:
                sec_monomer_group = _create_monomer_group(
                    restype, sec_molecule, particle_arrays
                )
                sec_molecule.sub_systems.append(sec_monomer_group)

        particle_indices = np.where(particle_arrays['molnums'] == molecule)[0]

        mol_system = self._create_system_node(
            name=molecule,
            branch_label='molecule',
            particle_indices=particle_indices,
        )

        # Check if molecule has multiple residues
        mol_resids = np.unique(particle_arrays['resids'][mol_system.particle_indices])
        if len(mol_resids) > 1:
            _add_residue_hierarchy(mol_system, particle_arrays)

        return mol_system

    def _parse_molecular_hierarchy(self, simulation: Simulation) -> None:
        """Parse molecular hierarchy (molecule groups, molecules, residues)."""

        def _get_particles_info() -> dict | None:
            """Get particle information from the first frame."""
            first_frame = 0
            particles_info = self._mdanalysistraj_parser.get('atoms_info', None)

            if particles_info is None:
                particles_info = self.traj_parsers.eval('atoms_info')
                if isinstance(particles_info, list):
                    particles_info = (
                        particles_info[first_frame]
                        if len(particles_info) > first_frame
                        else None
                    )

            return particles_info

        def _extract_particle_arrays(
            particles_info: dict, simulation: Simulation
        ) -> dict:
            """Extract and process particle information arrays."""
            first_frame = 0

            particle_labels = [
                ps.label for ps in simulation.model_system[first_frame].particle_states
            ]

            particles_elements = np.array(
                particles_info.get('elements', ['CGX'] * self.n_particles)
            )
            particles_types = np.array(particles_info.get('types', []))

            # Replace CGX placeholder elements if better labels available
            if 'CGX' in particles_elements:
                if particle_labels and 'CGX' not in particle_labels:
                    particles_elements = np.array(particle_labels)
                else:
                    particles_elements = particles_types

            return {
                'moltypes': np.array(particles_info.get('moltypes', [])),
                'molnums': np.array(particles_info.get('molnums', [])),
                'resids': np.array(particles_info.get('resids', [])),
                'resnames': np.array(particles_info.get('resnames', [])),
                'elements': particles_elements,
                'types': particles_types,
            }

        def _create_molecule_group(moltype: str, particle_arrays: dict) -> ModelSystem:
            """Create a molecule group with its constituent molecules."""
            particle_indices = np.where(particle_arrays['moltypes'] == moltype)[0]

            # Calculate composition formula
            mol_nums = particle_arrays['molnums'][particle_indices]
            moltype_count = np.unique(mol_nums).shape[0]

            molecule_group = self._create_system_node(
                name=f'group_{moltype}',
                branch_label='molecule_group',
                particle_indices=particle_indices,
                composition_formula=f'{moltype}({moltype_count})',
            )

            # Add individual molecules
            molecules = particle_arrays['molnums']
            for i_molecule, molecule in enumerate(
                np.unique(molecules[molecule_group.particle_indices])
            ):
                mol = self._create_molecule(molecule, i_molecule, particle_arrays)
                molecule_group.sub_systems.append(mol)

            return molecule_group

        particles_info = _get_particles_info()
        if particles_info is None:
            return
        particle_arrays = _extract_particle_arrays(particles_info, simulation)

        # Build molecular hierarchy
        moltypes = np.unique(particle_arrays['moltypes'])
        for moltype in moltypes:
            molecule_group = _create_molecule_group(moltype, particle_arrays)
            simulation.model_system[0].sub_systems.append(molecule_group)

    def parse_system(self, simulation):
        """Parse system information from trajectory and create model systems."""
        # Validate and prepare trajectory data
        if not self._validate_trajectory_data():
            return

        # Parse trajectory frames (dimension, positions, velocities, cell)
        self._parse_trajectory_frames(simulation)

        # Parse molecular hierarchy (molecule groups, molecules, residues)
        self._parse_molecular_hierarchy(simulation)

        # Mark the last (minimized/equilibrated) configuration as is_representative
        if simulation.model_system:
            simulation.model_system[-1].is_representative = True

    def _parse_input_structure(self, simulation: Simulation) -> None:
        """Store a known input structure when no trajectory was found.

        This is the structure supplied to LAMMPS by ``read_data``, not an
        inferred post-run or post-minimization configuration.
        """
        data_file = self._data_parser.mainfile
        if data_file is None or not os.path.isfile(data_file):
            self.logger.warning('No readable trajectory or input data file found.')
            return

        read_data_commands = []
        for log_parser in (self._log_parser, self._aux_log_parser):
            commands = log_parser.get('read_data') or []
            read_data_commands.extend(
                [commands] if isinstance(commands, str) else commands
            )
        if not any(command != 'CPU' for command in read_data_commands):
            self.logger.warning(
                'An input data file was found without a matching read_data command.'
            )
            return

        atom_styles = self._log_parser.get('atom_style') or self._aux_log_parser.get(
            'atom_style'
        )
        if not atom_styles or atom_styles[-1] != 'atomic':
            self.logger.warning(
                'Cannot reconstruct an input structure for LAMMPS atom style '
                f'{atom_styles}.'
            )
            return

        data_reader = MDAnalysisParser(format='DATA', atom_style='id type x y z')
        data_reader.logger = self.logger
        data_reader.mainfile = data_file
        try:
            data_reader.parse()
            frame = data_reader.get_frame_data(0)
            if frame['positions'] is None:
                self.logger.warning('Input data file has no readable atom positions.')
                return

            # MDAnalysis assigns angstrom to DATA coordinates regardless of the
            # LAMMPS unit style. The numeric values have the log's distance unit.
            distance_unit = self._log_parser.units.get('distance')
            if distance_unit is None:
                self.logger.warning('Cannot resolve the input structure distance unit.')
                return
            lattice_vectors = frame['lattice_vectors']
            if lattice_vectors is not None:
                lattice_vectors = lattice_vectors.magnitude * distance_unit

            boundary = self._log_parser.get('boundary') or self._aux_log_parser.get(
                'boundary'
            )
            boundary = boundary[-1] if boundary else ['p', 'p', 'p']
            particle_types = data_reader.get('atoms_info', {}).get('types')
            self.parse_trajectory_step(
                {
                    'lattice_vectors': lattice_vectors,
                    'periodic_boundary_conditions': [
                        value == 'p' for value in boundary
                    ],
                    # Type identifiers are labels, not evidence of elements.
                    'labels': particle_types,
                    'n_particles': data_reader.get_n_atoms(0),
                    'positions': frame['positions'].magnitude * distance_unit,
                    'dimensions': 3,
                },
                simulation,
                model_system=ModelSystem(
                    name='Input structure', is_representative=False
                ),
            )
            self.logger.warning(
                'Stored the LAMMPS input structure; no final structure is available.'
            )
        finally:
            data_reader.close()

    def parse_input(self, simulation: Simulation) -> None:
        """
        Parse input/control parameters from log file.

        TODO: Migrate from legacy parser
        - Extract input/output file information:
          * Data file basename
          * Trajectory file basename
        - Parse control parameters from log file commands
        - Map LAMMPS commands to x_lammps_inout_control_* attributes
        - Store in custom section x_lammps_section_control_parameters

        Legacy implementation: lines 1625-1650 in
        atomisticparsers/lammps/parser.py
        Key changes needed:
        - Determine if custom x_lammps sections should be migrated or dropped
        - Consider storing control parameters in standard schema if applicable
        - Update file path handling to use self._data_parser, self.traj_parsers
        """
        pass

    def _thermo_systems_by_step(
        self, simulation: Simulation
    ) -> dict[int, tuple[int, ModelSystem]]:
        """Match sampled frames by source step, excluding ambiguous duplicates."""
        if len(self.traj_parsers._parsers) != 1:
            return {}
        system_by_step: dict[int, tuple[int, ModelSystem]] = {}
        ambiguous_steps: set[int] = set()
        frame_indices = range(
            0, self._trajectory_frame_count, self.archive_sampling_rate
        )
        if len(frame_indices) == len(simulation.model_system):
            for frame_index, system in zip(
                frame_indices, simulation.model_system, strict=True
            ):
                step = system.time_step
                if step is None:
                    continue
                if step in system_by_step:
                    ambiguous_steps.add(step)
                else:
                    system_by_step[step] = (frame_index, system)
        for step in ambiguous_steps:
            system_by_step.pop(step, None)
        return system_by_step

    @staticmethod
    def _set_thermo_energies(
        output: WorkflowOutputs,
        values: dict[str, tuple[str, float]],
        energy_unit: Any,
    ) -> None:
        total = values.get('toteng')
        if total is None:
            if 'poteng' in values:
                output.potential_energies.append(
                    PotentialEnergy(
                        name=values['poteng'][0],
                        value=values['poteng'][1] * energy_unit,
                    )
                )
            if 'kineng' in values:
                output.kinetic_energies.append(
                    KineticEnergy(
                        name=values['kineng'][0],
                        value=values['kineng'][1] * energy_unit,
                    )
                )
            for key, (name, value) in values.items():
                if key.startswith('e_'):
                    output.potential_energies.append(
                        PotentialEnergy(name=name, value=value * energy_unit)
                    )
            return

        energy = TotalEnergy(value=total[1] * energy_unit)
        for key, contribution_type in (
            ('poteng', 'potential'),
            ('kineng', 'kinetic'),
        ):
            if key in values:
                name, value = values[key]
                energy.contributions.append(
                    BaseEnergy(
                        name=name,
                        contribution_type=contribution_type,
                        value=value * energy_unit,
                    )
                )
        for key, (name, value) in values.items():
            if key.startswith('e_'):
                contribution_type = _ENERGY_CONTRIBUTION_TYPES.get(key, key[2:])
                energy.contributions.append(
                    BaseEnergy(
                        name=name,
                        contribution_type=contribution_type,
                        value=value * energy_unit,
                    )
                )
        output.total_energies.append(energy)

    @staticmethod
    def _set_thermo_properties(
        output: WorkflowOutputs,
        values: dict[str, tuple[str, float]],
        units: dict[str, Any],
    ) -> set[str]:
        """Map typed scalar thermo values and report those without a unit."""
        unsupported: set[str] = set()
        mappings = (
            ('press', 'pressures', Pressure, units.get('pressure')),
            ('volume', 'volumes', Volume, units['distance'] ** 3),
            ('density', 'mass_densities', MassDensity, units.get('density')),
            ('enthalpy', 'enthalpies', Enthalpy, units['energy']),
        )
        for key, section_name, section_type, unit in mappings:
            if key not in values:
                continue
            name, value = values[key]
            if unit is None:
                unsupported.add(name)
                continue
            getattr(output, section_name).append(section_type(value=value * unit))
        return unsupported

    def _set_thermo_forces(
        self,
        output: WorkflowOutputs,
        frame_index: int,
        system: ModelSystem,
    ) -> None:
        """Store complete force arrays only for an unambiguous source frame."""
        raw_forces = self.traj_parsers.eval('get_forces', frame_index)
        if raw_forces is None:
            return
        forces = self.apply_unit(raw_forces, 'force')
        values = np.asarray(forces.magnitude)
        if values.shape != (system.n_particles, 3) or not np.isfinite(values).all():
            self.logger.warning(
                'Skipping LAMMPS force array with inconsistent particle shape.'
            )
            return
        output.total_forces.append(TotalForce(value=forces))

    def parse_thermodynamic_data(  # noqa: PLR0912, PLR0915
        self, simulation: Simulation
    ) -> None:
        """Store the last thermo block as ordered, source-backed outputs.

        The current workflow parser selects one final run. Earlier thermo
        blocks are not folded into that run until multi-run workflows are
        supported. Unknown columns, including pressure, remain explicit gaps.
        """
        source = None
        blocks = []
        for reader in (self._aux_log_parser, self._log_parser):
            blocks = reader.get_thermodynamic_blocks()
            if blocks:
                source = reader
                break
        if source is None:
            return
        if len(blocks) > 1:
            self.logger.warning(
                'Multiple LAMMPS thermo blocks found; only the last is stored.',
                n_blocks=len(blocks),
            )
        records = blocks[-1]
        units = source.units
        system_by_step = self._thermo_systems_by_step(simulation)

        step_counts: dict[int, int] = {}
        for record in records:
            step = record.get('Step')
            if step is not None and np.isfinite(step) and float(step).is_integer():
                step = int(step)
                step_counts[step] = step_counts.get(step, 0) + 1

        steps = [int(record['Step']) for record in records if 'Step' in record]
        has_known_origin = (
            bool(steps)
            and steps[0] == 0
            and all(earlier < later for earlier, later in zip(steps, steps[1:]))
        )
        timestep = (
            self.archive.workflow2.method.integration_timestep
            if isinstance(self.archive.workflow2, MolecularDynamics)
            else None
        )
        output_type = (
            TrajectoryOutputs
            if isinstance(self.archive.workflow2, MolecularDynamics)
            else WorkflowOutputs
        )
        energy_unit = units['energy']
        temperature_unit = units['temperature']
        unsupported_columns: set[str] = set()
        max_step = np.iinfo(np.int32).max

        for record in records:
            raw_step = record.get('Step')
            if (
                raw_step is None
                or not np.isfinite(raw_step)
                or not float(raw_step).is_integer()
                or raw_step < 0
                or raw_step > max_step
            ):
                self.logger.warning('Skipping LAMMPS thermo row with invalid step.')
                continue
            step = int(raw_step)
            values = {key.lower(): (key, value) for key, value in record.items()}
            output = output_type(step=step)
            if step_counts[step] == 1 and step in system_by_step:
                frame_index, system = system_by_step[step]
                output.model_system_ref = system
                self._set_thermo_forces(output, frame_index, system)
            if isinstance(output, TrajectoryOutputs):
                if 'time' in values:
                    output.time = values['time'][1] * units['time']
                elif has_known_origin and timestep is not None:
                    output.time = step * timestep

            if 'temp' in values:
                output.temperatures.append(
                    Temperature(value=values['temp'][1] * temperature_unit)
                )
            if 'cpu' in values:
                # LAMMPS CPU is elapsed seconds since the start of this run,
                # independent of the selected physical unit style.
                output.cpu1_end = values['cpu'][1] * ureg.second

            self._set_thermo_energies(output, values, energy_unit)
            unsupported_columns.update(
                self._set_thermo_properties(output, values, units)
            )

            supported = {
                'step',
                'temp',
                'toteng',
                'poteng',
                'kineng',
                'cpu',
                'press',
                'volume',
                'density',
                'enthalpy',
            }
            if isinstance(output, TrajectoryOutputs):
                supported.add('time')
            unsupported_columns.update(
                name
                for key, (name, _) in values.items()
                if key not in supported and not key.startswith('e_')
            )
            simulation.outputs.append(output)

        if unsupported_columns:
            self.logger.warning(
                'LAMMPS thermo columns have no mapped output destination.',
                columns=sorted(unsupported_columns),
            )

    def parse_workflow(self, simulation: Simulation) -> None:
        """Attach method settings for the last MD run echoed in the log files."""

        def log_lines():
            for log_parser in (self._log_parser, self._aux_log_parser):
                if log_parser.mainfile and os.path.isfile(log_parser.mainfile):
                    with open(
                        log_parser.mainfile, encoding='utf-8', errors='replace'
                    ) as source:
                        yield from source

        settings = selected_run(log_lines())
        if settings is None or settings['kind'] != 'md':
            return

        # The preflight resolved and validated units from this log, including
        # runtime-only "Unit style" records when no command is echoed.
        if self._resolved_unit_style is None:
            raise RuntimeError('LAMMPS unit style must be resolved before parsing.')
        self.archive.workflow2 = MolecularDynamics(
            method=build_md_method(settings, get_unit(self._resolved_unit_style)),
        )

    def _configure_parsers(self) -> None:
        """Configure all parsers with loggers and basic settings."""
        # Configure main log parser
        self._log_parser.mainfile = self.mainfile
        self._log_parser.logger = self.logger
        if self._resolved_unit_style is None:
            raise RuntimeError('LAMMPS unit style must be resolved before setup.')
        self._log_parser._units = get_unit(self._resolved_unit_style)

        # Set up auxiliary log parser if specified
        aux_log_files = self._log_parser.get('log')
        if aux_log_files:
            self._aux_log_parser.mainfile = os.path.join(
                self._log_parser.maindir,
                aux_log_files[0],
            )
            # We assign units here which is read from log parser
            self._aux_log_parser._units = self._log_parser.units
            self._aux_log_parser.logger = self.logger

        # Configure trajectory parsers
        self._traj_parser.logger = self.logger
        self._traj_parser._chemical_symbols = None
        self._xyztraj_parser.logger = self.logger
        self._mdanalysistraj_parser.logger = self.logger

        # Configure data parser
        self._data_parser.logger = self.logger

    def _set_data_files(self) -> None:
        """Parse and configure data file(s) associated with calculation."""
        log_parser = self._log_parser
        if (
            not log_parser.get('read_data')
            and self._aux_log_parser.mainfile
            and os.path.isfile(self._aux_log_parser.mainfile)
            and self._aux_log_parser.get('read_data')
        ):
            log_parser = self._aux_log_parser
        data_files = log_parser.get_data_files()

        if len(data_files) > 1:
            self.logger.warning('Multiple data files are specified')

        if data_files:
            self._data_parser.mainfile = data_files[0]

    def _create_trajectory_parser(
        self,
        traj_file: str,
        index: int,
        data_file: str,
        dump_commands: list | None = None,
    ) -> TrajParser | XYZTrajParser | MDAnalysisParser:
        """
        Create appropriate trajectory parser based on file type.

        Parser initialization for each traj file cannot be avoided as there are
        cases where traj files can share the same parser.
        """

        def _create_formatted_parser(
            traj_file: str, file_type: str, data_file: str
        ) -> MDAnalysisParser:
            """Create MDAnalysis parser for specified trajectory file formats."""
            traj_parser = MDAnalysisParser(
                topology_format='DATA', format=file_type.upper()
            )
            traj_parser.mainfile = data_file
            traj_parser.auxilliary_files = [traj_file]
            self._mdanalysistraj_parser = traj_parser
            return traj_parser

        # TODO: Handling of file_type = 'atom' is a LB edit, test
        def _create_custom_parser(
            traj_file: str, index: int, data_file: str, file_type: str
        ) -> TrajParser | MDAnalysisParser:
            """Create parser for custom or atom LAMMPS dump formats."""
            custom_options = None
            if file_type == 'custom':
                custom_options = dump_commands[index][5:]
                # Convert unwrapped coordinates (xu, yu, zu) to regular (x, y, z)
                custom_options = [
                    option.replace('xu', 'x').replace('yu', 'y').replace('zu', 'z')
                    for option in custom_options
                ]
                custom_options = ' '.join(custom_options)

            # Try MDAnalysis first
            traj_parser = MDAnalysisParser(
                topology_format='DATA',
                format='LAMMPSDUMP',
                atom_style=custom_options,
            )
            traj_parser.mainfile = data_file
            traj_parser.auxilliary_files = [traj_file]

            # Check if MDAnalysis can construct the universe or parse the atoms,
            # otherwise will fall back to TrajParser
            if traj_parser.universe is None or 'CGX' in traj_parser.get(
                'atoms_info', {}
            ).get('names', []):
                # MDAnalysis is necessary to calculate rdf and atomsgroup
                if index == 0:
                    self._mdanalysistraj_parser = traj_parser
                traj_parser = TrajParser()
                traj_parser.mainfile = traj_file

            return traj_parser

        # Determine file type from dump command or file extension
        if dump_commands is None:
            dump_commands = self._log_parser.get('dump')
        if dump_commands:
            file_type = dump_commands[index][2]
        else:
            # TODO: Assumes the extension is always a valid lammps dump format, improve
            # Fallback to file extension
            file_type = traj_file.rsplit('.', maxsplit=1)[-1]

        # TODO: add support for other LAMMPS dump file formats (https://docs.lammps.org/dump.html)
        if file_type in ('dcd', 'xyz') and data_file:
            return _create_formatted_parser(traj_file, file_type, data_file)

        # TODO: 'atom' keyword is a LB edit, test
        elif file_type in ('custom', 'atom') and data_file:
            return _create_custom_parser(traj_file, index, data_file, file_type)

        else:
            self.logger.warning('File type of %s not recognized.', traj_file)
            traj_parser = TrajParser()
            traj_parser.mainfile = traj_file
            # TODO: provide support for other file types
            return traj_parser

    def _parse_trajectory_files(
        self,
    ) -> list[TrajParser | XYZTrajParser | MDAnalysisParser]:
        """Parse trajectory files and create appropriate parsers."""
        log_parser = self._log_parser
        if (
            not log_parser.get('dump')
            and self._aux_log_parser.mainfile
            and os.path.isfile(self._aux_log_parser.mainfile)
            and self._aux_log_parser.get('dump')
        ):
            log_parser = self._aux_log_parser
        traj_files = log_parser.get_traj_files()
        dump_commands = log_parser.get('dump')

        if len(traj_files) > 1:
            self.logger.warning('Multiple traj files are specified')

        data_file = self._data_parser.mainfile
        parsers = []

        for n, traj_file in enumerate(traj_files):
            if not os.path.isfile(traj_file):
                self.logger.warning(
                    f'Referenced trajectory file is missing: {traj_file}'
                )
                continue
            traj_parser = self._create_trajectory_parser(
                traj_file, n, data_file, dump_commands
            )
            parsers.append(traj_parser)

        self.traj_parsers = TrajParsers(parsers)
        return parsers

    def _parse_content_sections(self) -> None:
        self.parse_method(self.archive.data)
        self.parse_system(self.archive.data)
        if self.archive.data.model_method and self.archive.data.model_system:
            method = self.archive.data.model_method[-1]
            if isinstance(method, ForceField):
                atoms_info = self.traj_parsers.eval('atoms_info')
                if isinstance(atoms_info, list):
                    atoms_info = atoms_info[0] if atoms_info else None
                atom_ids = (
                    atoms_info.get('id') if isinstance(atoms_info, dict) else None
                )
                attach_particle_parameters(
                    method,
                    self._log_parser,
                    self._data_parser,
                    atom_ids,
                    self.archive.data.model_system[0].particle_states,
                    get_unit(self._resolved_unit_style),
                )

        # TODO: uncomment when implemented
        # self.parse_input(self.archive.data)

    def write_to_archive(self) -> None:
        self._bond_list = None
        self._trajectory_frame_count = 0
        self._resolved_unit_style = self._require_physical_unit_style()
        self.archive.data = Simulation(program=Program(name='LAMMPS'))
        # LAMMPS mainfile is the main log file
        self.basename = os.path.basename(self.mainfile)
        self.basedir = os.path.dirname(self.mainfile)

        # Configure all parsers (loggers, units, etc.)
        self._configure_parsers()

        # Set up and parse data files
        self._set_data_files()

        # Set up and parse trajectory files
        parsers = self._parse_trajectory_files()
        self.parse_workflow(self.archive.data)

        # Parse system, method, parameters, thermodynamic data, etc.
        if parsers:
            self._parse_content_sections()
        else:
            self._parse_input_structure(self.archive.data)

        self.parse_thermodynamic_data(self.archive.data)

        # Close all parser instances
        self._mdanalysistraj_parser.close()
        self._log_parser.close()
        self._aux_log_parser.close()
        self._data_parser.close()
        for parser in parsers:
            parser.close()


class LammpsParser(MatchingParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.archive_writer = LammpsArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
