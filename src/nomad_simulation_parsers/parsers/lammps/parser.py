import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    BarostatParameters,
    MolecularDynamicsMethod,
    ThermostatParameters,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.lammps.file_parsers import DataParser, LogParser
from nomad_simulation_parsers.parsers.lammps.trajectory_parsers import (
    TrajParser,
    TrajParsers,
    XYZTrajParser,
)
from nomad_simulation_parsers.parsers.utils.mdanalysisparser import MDAnalysisParser
from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser


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
        self._bond_list = None

    _magic_three: int = 3
    _magic_four: int = 4
    _magic_five: int = 5

    def apply_unit(self, value: Any, unit: str) -> float:
        if not hasattr(value, 'units'):
            value = value * self._log_parser.units.get(unit, 1)
        return value

    def _extract_temp_params(
        self, fix_cmd: list[str], temp_idx: int, coupling_idx: int
    ) -> tuple[Any, Any] | None:
        """Extract temperature and coupling constant from fix command."""
        if temp_idx + coupling_idx >= len(fix_cmd):
            return None
        try:
            temp = self.apply_unit(float(fix_cmd[temp_idx]), 'temperature')
            coupling = self.apply_unit(float(fix_cmd[coupling_idx]), 'time')
            _temp_params = (temp, coupling)
        except (ValueError, IndexError, TypeError):
            _temp_params = None
        return _temp_params

    def _extract_thermostat_settings(
        self, fix_commands: list[list[str]] | None
    ) -> ThermostatParameters | None:
        """
        Extract thermostat parameters from LAMMPS fix commands.

        Supports fix nvt, fix langevin, fix temp/berendsen, fix temp/rescale.

        Args:
            fix_commands: List of fix command arguments from log parser

        Returns:
            ThermostatParameters object or None if no thermostat found
        """
        if fix_commands is None:
            return None
        thermostat_map = {
            'nvt': 'nose_hoover',
            'npt': 'nose_hoover',
            'langevin': 'langevin_leap_frog',
            'temp/berendsen': 'berendsen',
            'temp/rescale': 'velocity_rescaling',
        }
        for fix_cmd in fix_commands:
            if len(fix_cmd) < self._magic_four:
                continue
            fix_type = fix_cmd[2]
            if fix_type not in thermostat_map:
                continue
            params = ThermostatParameters()
            params.thermostat_type = thermostat_map[fix_type]
            if fix_type in ['nvt', 'npt'] and 'temp' in fix_cmd:
                temp_idx = fix_cmd.index('temp')
                temp_params = self._extract_temp_params(
                    fix_cmd, temp_idx + 1, temp_idx + 3
                )
                if temp_params:
                    params.reference_temperature, params.coupling_constant = temp_params
            elif fix_type in ['langevin', 'temp/berendsen']:
                temp_params = self._extract_temp_params(fix_cmd, 3, 5)
                if temp_params:
                    params.reference_temperature, params.coupling_constant = temp_params
            elif fix_type == 'temp/rescale' and len(fix_cmd) > self._magic_four:
                try:
                    params.reference_temperature = self.apply_unit(
                        float(fix_cmd[4]), 'temperature'
                    )
                except (ValueError, TypeError):
                    pass
            return params
        return None

    def _extract_pressure_params(
        self, fix_cmd: list[str], style_idx: int, offset: int
    ) -> tuple[Any, Any] | None:
        """Extract pressure and coupling matrices from fix command."""
        if style_idx + offset >= len(fix_cmd):
            return None
        try:
            press = self.apply_unit(float(fix_cmd[style_idx + 1]), 'pressure')
            coupling = self.apply_unit(float(fix_cmd[style_idx + offset]), 'time')
            _pressure_params = (np.eye(3) * press, np.eye(3) * coupling)
        except (ValueError, IndexError, TypeError):
            _pressure_params = None
        return _pressure_params

    def _extract_modulus(self, fix_cmd: list[str]) -> np.ndarray | None:
        """Extract modulus (bulk modulus) from fix command and convert to
        compressibility.

        LAMMPS uses 'modulus' keyword in fix npt/nph commands to specify bulk modulus.
        Compressibility = 1 / modulus (in consistent units).

        Args:
            fix_cmd: Fix command arguments list

        Returns:
            3x3 compressibility matrix (1/pascal) or None if modulus not found
        """
        result = None
        if 'modulus' not in fix_cmd:
            return result

        try:
            modulus_idx = fix_cmd.index('modulus')
            # Modulus value follows the 'modulus' keyword
            modulus_value = float(fix_cmd[modulus_idx + 1])

            # Apply pressure units to modulus (same units as pressure in LAMMPS)
            modulus_with_units = self.apply_unit(modulus_value, 'pressure')

            # Compressibility = 1 / modulus
            # Units: 1/pascal = 1 / (kg/(m·s²)) = m·s²/kg = Pa⁻¹
            compressibility_value = 1.0 / modulus_with_units.magnitude
            compressibility_unit = 1.0 / modulus_with_units.units

            # Create 3x3 diagonal matrix (isotropic compressibility)
            result = np.eye(3) * compressibility_value * compressibility_unit
        except (ValueError, IndexError, TypeError, ZeroDivisionError):
            result = None

        return result

    def _extract_timestep(self) -> Any | None:
        """Extract integration timestep from log parser."""
        timestep_value = self._log_parser.get('timestep')
        if timestep_value is None:
            return None
        if isinstance(timestep_value, list):
            timestep_value = timestep_value[-1]
        try:
            _timestep = self.apply_unit(float(timestep_value), 'time')
        except (ValueError, TypeError):
            _timestep = None
        return _timestep

    def _extract_n_steps(self) -> int | None:
        """Extract number of steps from run command."""
        run_command = self._log_parser.get('run')
        if run_command is None:
            return None
        if isinstance(run_command, list) and isinstance(run_command[0], list):
            run_command = run_command[-1]
        try:
            _n_steps = int(
                run_command[0] if isinstance(run_command, list) else run_command
            )
        except (ValueError, TypeError, IndexError):
            _n_steps = None
        return _n_steps

    def _extract_thermo_frequency(self) -> int | None:
        """Extract thermodynamics output frequency."""
        thermo_freq = self._log_parser.get('thermo')
        if thermo_freq is None:
            return None
        if isinstance(thermo_freq, list):
            thermo_freq = thermo_freq[-1]
        try:
            _thermo_freq = int(thermo_freq)
        except (ValueError, TypeError):
            _thermo_freq = None
        return _thermo_freq

    def _extract_masses(self) -> np.ndarray | None:
        """Extract particle masses from data file."""
        masses_data = self._data_parser.get('Masses', None)
        if (
            not masses_data
            or not isinstance(masses_data, list)
            or len(masses_data) == 0
        ):
            return None
        return masses_data[0][1] if isinstance(masses_data[0], tuple) else masses_data

    def _set_parser_masses(self, masses: np.ndarray | None) -> None:
        """Set masses on trajectory parsers."""
        for parser in self.traj_parsers._parsers:
            if isinstance(parser, TrajParser):
                parser.masses = masses

    def _extract_integrator_type(self) -> str | None:
        """Extract integrator type from run_style command.

        LAMMPS run_style options: verlet, verlet/split, respa, respa/omp
        Maps to NOMAD integrator_type enum.

        Returns:
            Integrator type string or None if not found
        """
        run_style = self._log_parser.get('run_style')
        if run_style is None:
            return None

        # Handle nested list structure from multiple run_style commands
        # E.g., [['verlet'], ['respa', '4', '2', ...]] -> take last command
        if isinstance(run_style, list) and len(run_style) > 0:
            if isinstance(run_style[0], list):
                # Nested list: take last command
                run_style = run_style[-1]

        # Extract style name from command arguments
        # E.g., ['respa', '4', '2', ...] -> 'respa'
        #       'verlet' -> 'verlet'
        if isinstance(run_style, list) and len(run_style) > 0:
            style = run_style[0]
        else:
            style = run_style

        if style is None:
            return None

        result = None
        style_lower = str(style).lower()

        # Map LAMMPS run_style to NOMAD integrator_type
        integrator_map = {
            'verlet': 'velocity_verlet',
            'verlet/split': 'velocity_verlet_split',
            'respa': 'respa',  # reversible reference system propagator algorithm
            'respa/omp': 'respa',
        }

        result = integrator_map.get(style_lower, 'velocity_verlet')  # default fallback
        return result

    def _extract_barostat_settings(
        self, fix_commands: list[list[str]] | None
    ) -> BarostatParameters | None:
        """
        Extract barostat parameters from LAMMPS fix commands.

        Supports fix npt, fix nph, fix press/berendsen.

        Args:
            fix_commands: List of fix command arguments from log parser

        Returns:
            BarostatParameters object or None if no barostat found
        """
        if fix_commands is None:
            return None
        coupling_styles = {
            'iso': 'isotropic',
            'aniso': 'anisotropic',
            'tri': 'anisotropic',
        }
        for fix_cmd in fix_commands:
            if len(fix_cmd) < self._magic_four:
                continue
            fix_type = fix_cmd[2]
            if fix_type not in ['npt', 'nph', 'press/berendsen']:
                continue
            params = BarostatParameters()
            params.barostat_type = (
                'nose_hoover' if fix_type in ['npt', 'nph'] else 'berendsen'
            )
            offset = 3 if fix_type in ['npt', 'nph'] else 2
            for style, coupling_type in coupling_styles.items():
                if style not in fix_cmd:
                    continue
                params.coupling_type = coupling_type
                style_idx = fix_cmd.index(style)
                pressure_params = self._extract_pressure_params(
                    fix_cmd, style_idx, offset
                )
                if pressure_params:
                    params.reference_pressure, params.coupling_constant = (
                        pressure_params
                    )
                break
            # Extract modulus and convert to compressibility
            compressibility = self._extract_modulus(fix_cmd)
            if compressibility is not None:
                params.compressibility = compressibility

            # Extract modulus and convert to compressibility
            compressibility = self._extract_modulus(fix_cmd)
            if compressibility is not None:
                params.compressibility = compressibility

            return params
        return None

    def _determine_ensemble(
        self, has_thermostat: bool, has_barostat: bool
    ) -> str | None:
        """
        Determine thermodynamic ensemble from presence of thermostat and barostat.

        Args:
            has_thermostat: Whether thermostat parameters were found
            has_barostat: Whether barostat parameters were found

        Returns:
            Ensemble type string ('NVE', 'NVT', 'NPH', 'NPT') or None
        """
        result = None
        if not has_thermostat and not has_barostat:
            result = 'NVE'
        elif has_thermostat and not has_barostat:
            result = 'NVT'
        elif not has_thermostat and has_barostat:
            result = 'NPH'
        elif has_thermostat and has_barostat:
            result = 'NPT'
        return result

    def _extract_dump_frequencies(
        self, dump_commands: list[list[str]] | None
    ) -> dict[str, int | None]:
        """
        Extract save frequencies from LAMMPS dump commands.

        Args:
            dump_commands: List of dump command arguments from log parser

        Returns:
            Dictionary with coordinate/velocity/force save frequencies
        """
        result = {
            'coordinate': None,
            'velocity': None,
            'force': None,
        }
        if dump_commands is None:
            return result
        for dump_cmd in dump_commands:
            if len(dump_cmd) < self._magic_five:
                continue
            try:
                frequency = int(dump_cmd[1])
            except (ValueError, IndexError):
                continue
            dump_vars = dump_cmd[5:] if len(dump_cmd) > self._magic_five else []
            has_coords = any(
                var in dump_vars for var in ['x', 'y', 'z', 'xu', 'yu', 'zu']
            )
            has_velocities = any(var in dump_vars for var in ['vx', 'vy', 'vz'])
            has_forces = any(var in dump_vars for var in ['fx', 'fy', 'fz'])
            if has_coords and result['coordinate'] is None:
                result['coordinate'] = frequency
            if has_velocities and result['velocity'] is None:
                result['velocity'] = frequency
            if has_forces and result['force'] is None:
                result['force'] = frequency
        return result

    def parse_method(self, simulation: Simulation) -> None:
        """
        Parse method information from data file and log file.

        Extracts:
        - Integration parameters (timestep, n_steps, integrator_type)
        - Thermostat parameters from fix commands (nvt, langevin, temp/berendsen, etc.)
        - Barostat parameters from fix commands (npt, nph, press/berendsen)
        - Thermodynamic ensemble (NVE, NVT, NPH, NPT)
        - Save frequencies from dump/thermo commands

        TODO: Still to migrate from legacy parser
        - Parse interactions (bonds, angles, dihedrals, impropers, pair_coeffs,
          bond_coeffs, angle_coeffs, etc.) using MDAnalysis
        - Set ForceField with Model containing interactions
        - Parse force calculation parameters:
          * pair_style: extract vdw_cutoff, coulomb_cutoff
          * kspace_style: set coulomb_type (ewald, particle_particle_particle_mesh,
            multilevel_summation)
        - Parse neighbor searching parameters:
          * neighbor: set neighbor_update_cutoff (add to vdw_cutoff)
          * neigh_modify: extract neighbor_update_frequency from 'every' parameter

        Legacy implementation: lines 1531-1624 in atomisticparsers/lammps/parser.py
        """

        if self.traj_parsers[0].mainfile is None or self._data_parser.mainfile is None:
            return

        if self.traj_parsers.eval('n_frames') is None:
            return

        method = MolecularDynamicsMethod()

        integrator_type = self._extract_integrator_type()
        if integrator_type is not None:
            method.integrator_type = integrator_type

        timestep = self._extract_timestep()
        if timestep is not None:
            method.integration_timestep = timestep

        n_steps = self._extract_n_steps()
        if n_steps is not None:
            method.n_steps = n_steps

        fix_commands = self._log_parser.get('fix')
        thermostat = self._extract_thermostat_settings(fix_commands)
        if thermostat is not None:
            method.thermostat_parameters = [thermostat]

        barostat = self._extract_barostat_settings(fix_commands)
        if barostat is not None:
            method.barostat_parameters = [barostat]

        ensemble = self._determine_ensemble(
            thermostat is not None, barostat is not None
        )
        if ensemble is not None:
            method.thermodynamic_ensemble = ensemble

        frequencies = self._extract_dump_frequencies(self._log_parser.get('dump'))
        if frequencies['coordinate'] is not None:
            method.coordinate_save_frequency = frequencies['coordinate']
        if frequencies['velocity'] is not None:
            method.velocity_save_frequency = frequencies['velocity']
        if frequencies['force'] is not None:
            method.force_save_frequency = frequencies['force']

        thermo_freq = self._extract_thermo_frequency()
        if thermo_freq is not None:
            method.thermodynamics_save_frequency = thermo_freq

        masses = self._extract_masses()
        self._set_parser_masses(masses)

        simulation.model_method.append(method)

    def _validate_trajectory_data(self) -> bool:
        """Validate that trajectory data is available and extract basic info."""
        n_traj = self.traj_parsers.eval('n_frames')
        if n_traj is None:
            return False

        self.n_atoms = [self.traj_parsers.eval('get_n_atoms', n) for n in range(n_traj)]
        self.trajectory_steps = [
            step
            for n in range(n_traj)
            if (step := self.traj_parsers.eval('get_step', n)) is not None
        ]
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
            else:
                # Convert from 1-based (LAMMPS data file) to 0-based indexing (NOMAD)
                # Extract columns 2:4 (atom IDs), convert to int, subtract 1
                bond_array = bonds[0][1][:, 2:4].astype(int) - 1
                self._bond_list = bond_array

        for step in self.trajectory_steps:
            traj_n = self.trajectory_steps.index(step)

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
                'labels': self.traj_parsers.eval('get_atom_labels', traj_n),
                'n_particles': self.traj_parsers.eval('get_n_atoms', traj_n),
                'positions': self.apply_unit(
                    self.traj_parsers.eval('get_positions', traj_n), 'distance'
                ),
                'velocities': velocities,
                'bond_list': self._bond_list if self._bond_list is not None else None,
                'dimensions': dimension,
            }
            self.parse_trajectory_step(particles_dict, simulation)

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
                particles_info.get('elements', ['CGX'] * self.n_atoms)
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

    def parse_thermodynamic_data(self, simulation: Simulation) -> None:
        """
        Parse thermodynamic output data from log file.

        TODO: Migrate from legacy parser
        - Extract thermodynamic data from log file or aux log file
        - Map thermodynamic quantities to TrajectoryOutputs:
          * Step number and physical time (step * timestep)
          * Energy contributions (kinetic, potential, pair, bond, angle, etc.)
          * Total energy (TotEng)
          * Pressure and temperature
          * Forces (if available in trajectory)
          * Calculation time (CPU)
        - Create TrajectoryOutputs for each thermodynamic step
        - Link outputs to corresponding model_system via model_system_ref

        Legacy implementation: lines 971-1036 in atomisticparsers/lammps/parser.py
        Key changes needed:
        - Use nomad_simulations TrajectoryOutputs instead of Calculation
        - Map energy types using self._energy_mapping
        - Coordinate with parse_output_step in mdparserutils.py (lines 186-238)
        - Extract timestep from log file (get_time_step method)
        - Match thermodynamics_steps with trajectory_steps
        """
        pass

    def parse_workflow(self, simulation: Simulation) -> None:
        """
        Parse workflow information for geometry optimization runs.

        TODO: Migrate from legacy parser
        - Detect minimization runs from log file (minimization_stats)
        - Create GeometryOptimization workflow:
          * Extract minimization method (cg, hftn, sd, quickmin, fire, spin)
          * Map to standard method names (polak_ribiere_conjugant_gradient,
            hessian_free_truncated_newton, steepest_descent, damped_dynamics)
          * Parse optimization parameters (max steps, force convergence)
        - Extract optimization results:
          * Final energy difference
          * Final force maximum
        - Parse minimize/minimize/kk parameters:
          * optimization_steps_maximum
          * convergence_tolerance_force_maximum
        - Handle unit conversions using self._log_parser.units

        Legacy implementation: lines 1037-1120 in
        atomisticparsers/lammps/parser.py
        Key changes needed:
        - Determine which simulation workflow schema to use from nomad_simulations
        - Check if GeometryOptimization is available in nomad_simulations
        - Adapt to simulation.workflow structure if different from sec_run.workflow
        - Update unit conversion to use apply_unit method
        """
        pass

    def _configure_parsers(self) -> None:
        """Configure all parsers with loggers and basic settings."""
        # Configure main log parser
        self._log_parser.mainfile = self.mainfile
        self._log_parser.logger = self.logger
        self._log_parser._units = None

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
        data_files = self._log_parser.get_data_files()

        if len(data_files) > 1:
            self.logger.warning('Multiple data files are specified')

        if data_files:
            self._data_parser.mainfile = data_files[0]

    def _create_trajectory_parser(
        self, traj_file: str, index: int, data_file: str
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
                custom_options = self._log_parser.get('dump')[index][5:]
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
        dump_commands = self._log_parser.get('dump')
        if dump_commands:
            file_type = dump_commands[index][2]
        else:
            # TODO: Assumes the extension is always a valid lammps dump format, improve
            # Fallback to file extension
            file_type = traj_file.split('.', 1)[-1]

        # TODO: add support for other LAMMPS dump file formats (https://docs.lammps.org/dump.html)
        if file_type == 'dcd' or file_type == 'xyz' and data_file:
            return _create_formatted_parser(traj_file, data_file, file_type)

        # TODO: 'atom' keyword is a LB edit, test
        elif file_type == 'custom' or file_type == 'atom' and data_file:
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
        traj_files = self._log_parser.get_traj_files()

        if len(traj_files) > 1:
            self.logger.warning('Multiple traj files are specified')

        data_file = self._data_parser.mainfile
        parsers = []

        for n, traj_file in enumerate(traj_files):
            traj_parser = self._create_trajectory_parser(traj_file, n, data_file)
            parsers.append(traj_parser)

        self.traj_parsers = TrajParsers(parsers)
        return parsers

    def _parse_content_sections(self) -> None:
        self.parse_method(self.archive.data)
        self.parse_system(self.archive.data)

        # TODO: uncomment when implemented
        # self.parse_input(self.archive.data)
        # self.parse_thermodynamic_data(self.archive.data)
        # self.parse_workflow(self.archive.data)

    def write_to_archive(self) -> None:
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
        if not self.traj_parsers or self.traj_parsers[0] is None:
            return

        # Parse system, method, parameters, thermodynamic data, etc.
        self._parse_content_sections()

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
