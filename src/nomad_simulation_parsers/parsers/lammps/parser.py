import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg as _ureg
from nomad_simulations.schema_packages.force_field import (
    ForceCalculations,
    ForceField,
    HarmonicAngle,
    HarmonicBond,
    Potential,
)
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    BarostatParameters,
    DiffusionConstant,
    MeanSquaredDisplacement,
    MolecularDynamics,
    MolecularDynamicsMethod,
    MolecularDynamicsResults,
    RadialDistributionFunction,
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

_KSPACE_COULOMB_TYPE = {
    'ewald': 'ewald',
    'pppm': 'particle_particle_particle_mesh',
    'msm': 'multilevel_summation',
    'pppm/stagger': 'particle_particle_particle_mesh',
    'pppm/cg': 'particle_particle_particle_mesh',
    'pppm/disp': 'particle_particle_particle_mesh',
}

_FIX_ENSEMBLE = {
    'nve': 'NVE',
    'nvt': 'NVT',
    'npt': 'NPT',
    'nph': 'NPH',
    'langevin': 'NVT',
}

_THERMOSTAT_TYPE = {
    'nvt': 'nose_hoover',
    'npt': 'nose_hoover',
    'langevin': 'langevin_schneider',
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

    def apply_unit(self, value: Any, unit: str) -> float:
        if not hasattr(value, 'units'):
            value = value * self._log_parser.units.get(unit, 1)
        return value

    def parse_method(self, simulation: Simulation) -> None:
        """Parse ForceField method with force calculation parameters."""
        if self.traj_parsers[0].mainfile is None or self._data_parser.mainfile is None:
            return

        if self.traj_parsers.eval('n_frames') is None:
            return

        masses_data = self._data_parser.get('Masses', None)
        if masses_data and isinstance(masses_data, list) and len(masses_data) > 0:
            masses = (
                masses_data[0][1] if isinstance(masses_data[0], tuple) else masses_data
            )
        else:
            masses = None

        for parser in self.traj_parsers._parsers:
            if isinstance(parser, TrajParser):
                parser.masses = masses

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

        force_field = ForceField()
        self._parse_force_calculations(force_field)
        self._parse_interactions(force_field)
        simulation.model_method.append(force_field)

    def _parse_force_calculations(self, force_field: ForceField) -> None:
        """Populate ForceCalculations from pair_style, kspace, and neighbor settings."""
        force_calcs = ForceCalculations()
        dist = self._log_parser.units.get('distance', 1)
        _set_pair_cutoffs(self._log_parser.get('pair_style'), force_calcs, dist)
        _set_kspace_coulomb(self._log_parser.get('kspace_style'), force_calcs)
        _set_neighbor_cutoff(self._log_parser.get('neighbor'), force_calcs, dist)
        _set_neighbor_frequency(self._log_parser.get('neigh_modify'), force_calcs)
        force_field.numerical_settings.append(force_calcs)

    def _parse_interactions(self, force_field: ForceField) -> None:
        """Add Potential contributions from per-category style getters."""
        dist = self._log_parser.units.get('distance', 1)
        energy = self._log_parser.units.get('energy', 1)

        bond_styles = self._log_parser.get('bond_style') or []
        bond_coeffs = self._log_parser.get('bond_coeff') or []
        for style in [bond_styles] if isinstance(bond_styles, str) else bond_styles:
            style_str = (
                style if isinstance(style, str) else ' '.join(str(s) for s in style)
            )
            if 'harmonic' in style_str.lower():
                force_field.contributions.append(
                    _build_harmonic_bond(style_str, bond_coeffs, dist, energy)
                )

        angle_styles = self._log_parser.get('angle_style') or []
        angle_coeffs = self._log_parser.get('angle_coeff') or []
        for style in [angle_styles] if isinstance(angle_styles, str) else angle_styles:
            style_str = (
                style if isinstance(style, str) else ' '.join(str(s) for s in style)
            )
            if 'harmonic' in style_str.lower():
                force_field.contributions.append(
                    _build_harmonic_angle(style_str, angle_coeffs, energy)
                )

        pair_styles = self._log_parser.get('pair_style') or []
        for entry in [pair_styles] if isinstance(pair_styles, str) else pair_styles:
            style_str = entry if isinstance(entry, str) else str(entry[0])
            style_lower = style_str.lower().split()[0]
            if 'lj' in style_lower:
                pot = Potential()
                pot.type = 'nonbonded'
                pot.functional_form = 'lennard-jones'
                pot.name = style_str
                force_field.contributions.append(pot)

    def _validate_trajectory_data(self) -> bool:
        """Validate that trajectory data is available and extract basic info."""
        n_traj = self.traj_parsers.eval('n_frames')
        if n_traj is None:
            return False

        self.n_particles = [
            self.traj_parsers.eval('get_n_atoms', n) for n in range(n_traj)
        ]
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
                # Convert from List[Tuple[None, np.ndarray]] to List[Tuple[int, int]]
                self._bond_list = list(map(tuple, bonds[0][1][:, 2:4].astype(int)))

        # dimension is frame-independent; read once
        dimension = self._log_parser.get('dimension', 3)

        for step in self.trajectory_steps:
            traj_n = self.trajectory_steps.index(step)
            velocities = _get_quantity_with_units(traj_n, 'velocities', 'velocity')

            if traj_n == 0:
                _extract_bond_list()
                particles_dict = {
                    'lattice_vectors': _get_quantity_with_units(
                        traj_n, 'lattice_vectors', 'distance'
                    ),
                    'periodic_boundary_conditions': self.traj_parsers.eval(
                        'get_pbc', traj_n
                    ),
                    'labels': self.traj_parsers.eval('get_atom_labels', traj_n),
                    'n_particles': self.traj_parsers.eval('get_n_atoms', traj_n),
                    'positions': self.apply_unit(
                        self.traj_parsers.eval('get_positions', traj_n), 'distance'
                    ),
                    'velocities': velocities,
                    'bond_list': self._bond_list if self._bond_list else None,
                    'dimensions': dimension,
                }
            else:
                particles_dict = {
                    'lattice_vectors': None,
                    'periodic_boundary_conditions': None,
                    'labels': None,
                    'dimensions': dimension,
                    'positions': self.apply_unit(
                        self.traj_parsers.eval('get_positions', traj_n), 'distance'
                    ),
                    'velocities': velocities,
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

    def parse_input(self, simulation: Simulation) -> None:
        pass

    def parse_thermodynamic_data(self, simulation: Simulation) -> None:
        """Parse thermodynamic outputs from the auxiliary (or main) log file."""
        thermo_data = self._aux_log_parser.get_thermodynamic_data()
        if thermo_data is None:
            thermo_data = self._log_parser.get_thermodynamic_data()
        if thermo_data is None:
            return

        steps_raw = thermo_data.get('Step')
        if steps_raw is None or len(steps_raw) == 0:
            return

        self.thermodynamics_steps = [int(s) for s in steps_raw]

        for i, step in enumerate(steps_raw):
            step_int = int(step)
            data: dict[str, Any] = {'step': step_int}

            if 'Temp' in thermo_data:
                data['temperatures'] = {'value': thermo_data['Temp'][i]}

            if 'TotEng' in thermo_data:
                data['total_energies'] = {'value': thermo_data['TotEng'][i]}

            if 'PotEng' in thermo_data:
                data['potential_energies'] = {'value': thermo_data['PotEng'][i]}

            if 'KinEng' in thermo_data:
                data['kinetic_energies'] = {'value': thermo_data['KinEng'][i]}

            self.parse_output_step(data, simulation)

    def parse_workflow(self, simulation: Simulation) -> None:
        """Parse MolecularDynamics workflow with method settings and RDF/MSD results."""
        sampling_method, fix_style = self._log_parser.get_sampling_method()
        if sampling_method != 'molecular_dynamics':
            return

        md_method = self._build_md_method(fix_style)
        md_results = self._build_md_results()

        sec_md = MolecularDynamics()
        sec_md.method = md_method
        sec_md.results = md_results
        self.archive.workflow2 = sec_md

    def _build_md_method(self, fix_style: str) -> MolecularDynamicsMethod:
        """Build MolecularDynamicsMethod from log parser settings."""
        md_method = MolecularDynamicsMethod()

        ensemble = _FIX_ENSEMBLE.get(fix_style.lower())
        if ensemble is not None:
            md_method.thermodynamic_ensemble = ensemble

        timestep_raw = self._log_parser.get('timestep')
        if timestep_raw is not None and len(timestep_raw) > 0:
            time_unit = self._log_parser.units.get('time', 1)
            md_method.integration_timestep = float(timestep_raw[0]) * time_unit

        thermo_data = self._aux_log_parser.get_thermodynamic_data()
        if thermo_data is None:
            thermo_data = self._log_parser.get_thermodynamic_data()
        if thermo_data is not None:
            steps_raw = thermo_data.get('Step')
            if steps_raw is not None and len(steps_raw) > 0:
                md_method.n_steps = int(steps_raw[-1])

        thermo_freq = self._log_parser.get('thermo')
        if thermo_freq is not None and len(thermo_freq) > 0:
            md_method.thermodynamics_save_frequency = int(thermo_freq[0])

        thermostat_settings = self._log_parser.get_thermostat_settings()
        self._add_thermostat_parameters(md_method, fix_style, thermostat_settings)
        self._add_barostat_parameters(md_method, fix_style, thermostat_settings)

        return md_method

    def _add_thermostat_parameters(
        self,
        md_method: MolecularDynamicsMethod,
        fix_style: str,
        thermostat_settings: dict,
    ) -> None:
        """Add ThermostatParameters if temperature coupling is present."""
        target_t = thermostat_settings.get('target_T')
        tau = thermostat_settings.get('thermostat_tau')
        if target_t is None and tau is None:
            return

        thermostat = ThermostatParameters()
        thermostat_type = _THERMOSTAT_TYPE.get(fix_style.lower())
        if thermostat_type is not None:
            thermostat.thermostat_type = thermostat_type
        if target_t is not None:
            thermostat.reference_temperature = target_t
        if tau is not None:
            thermostat.coupling_constant = tau
        md_method.thermostat_parameters.append(thermostat)

    def _add_barostat_parameters(
        self,
        md_method: MolecularDynamicsMethod,
        fix_style: str,
        thermostat_settings: dict,
    ) -> None:
        """Add BarostatParameters if pressure coupling is present."""
        target_p = thermostat_settings.get('target_P')
        tau_p = thermostat_settings.get('barostat_tau')
        if target_p is None and tau_p is None:
            return

        barostat = BarostatParameters()
        if target_p is not None:
            pressure_pa = target_p.to('pascal').magnitude
            ref_p = np.zeros((3, 3))
            ref_p[0, 0] = ref_p[1, 1] = ref_p[2, 2] = pressure_pa
            barostat.reference_pressure = ref_p
        if tau_p is not None:
            tau_s = tau_p.to('second').magnitude
            coupling = np.zeros((3, 3))
            coupling[0, 0] = coupling[1, 1] = coupling[2, 2] = tau_s
            barostat.coupling_constant = coupling
        md_method.barostat_parameters.append(barostat)

    def _build_md_results(self) -> MolecularDynamicsResults:
        """Build MolecularDynamicsResults with RDF and MSD from MDAnalysis."""
        md_results = MolecularDynamicsResults()

        if self._mdanalysistraj_parser.universe is None:
            return md_results

        rdf_data = self._mdanalysistraj_parser.calc_molecular_rdf()
        if rdf_data is not None:
            self._populate_rdf(md_results, rdf_data)

        msd_data = (
            self._mdanalysistraj_parser.calc_molecular_mean_squared_displacements()
        )
        if msd_data is not None:
            self._populate_msd(md_results, msd_data)

        return md_results

    def _populate_rdf(
        self, md_results: MolecularDynamicsResults, rdf_data: dict
    ) -> None:
        """Add RadialDistributionFunction entries to MolecularDynamicsResults."""
        types = rdf_data.get('types', [])
        bins_list = rdf_data.get('bins', [])
        value_list = rdf_data.get('value', [])
        frame_start_list = rdf_data.get('frame_start', [])
        frame_end_list = rdf_data.get('frame_end', [])
        n_smooth = rdf_data.get('n_smooth', 0)

        for i, pair_type in enumerate(types):
            if i >= len(bins_list) or i >= len(value_list):
                continue
            rdf = RadialDistributionFunction()
            rdf.label = str(pair_type)
            rdf.type = 'molecular'
            rdf.n_smooth = n_smooth
            rdf.bins = bins_list[i]
            rdf.value = value_list[i]
            rdf.n_bins = len(bins_list[i]) if bins_list[i] is not None else 0
            if i < len(frame_start_list):
                rdf.frame_start = int(frame_start_list[i])
            if i < len(frame_end_list):
                rdf.frame_end = int(frame_end_list[i])
            md_results.radial_distribution_functions.append(rdf)

    def _populate_msd(
        self, md_results: MolecularDynamicsResults, msd_data: dict
    ) -> None:
        """Add MeanSquaredDisplacement and DiffusionConstant entries."""
        types = msd_data.get('types', [])
        times_arr = msd_data.get('times')
        value_arr = msd_data.get('value')
        diffusion_arr = msd_data.get('diffusion_constant', [])
        error_arr = msd_data.get('error_diffusion_constant', [])

        for i, mol_type in enumerate(types):
            msd = MeanSquaredDisplacement()
            msd.label = str(mol_type)
            msd.type = 'molecular'
            msd.direction = 'xyz'
            if times_arr is not None and i < len(times_arr):
                msd.times = times_arr[i]
                msd.n_times = len(times_arr[i])
            if value_arr is not None and i < len(value_arr):
                msd.value = value_arr[i]
            md_results.mean_squared_displacements.append(msd)

            if i < len(diffusion_arr):
                dc = DiffusionConstant()
                dc.label = str(mol_type)
                dc.value = diffusion_arr[i]
                if i < len(error_arr):
                    dc.n_smooth = int(error_arr[i]) if error_arr[i] is not None else 0
                md_results.diffusion_constants.append(dc)

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
            file_type = traj_file.rsplit('.', maxsplit=1)[-1]

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
        self.parse_thermodynamic_data(self.archive.data)
        self.parse_workflow(self.archive.data)

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


def _set_pair_cutoffs(
    pair_styles: list | None,
    force_calcs: ForceCalculations,
    dist_unit: Any,
) -> None:
    """Set vdw_cutoff and coulomb_cutoff on force_calcs from pair_style data."""
    if not pair_styles:
        return
    pair_entry = pair_styles[0]
    _MIN_PAIR_ENTRY = 2
    if not (isinstance(pair_entry, list) and len(pair_entry) >= _MIN_PAIR_ENTRY):
        return
    style_name = str(pair_entry[0])
    try:
        cutoff = float(pair_entry[-1]) * dist_unit
        force_calcs.vdw_cutoff = cutoff
        if 'coul' in style_name:
            force_calcs.coulomb_cutoff = cutoff
    except (ValueError, TypeError):
        pass


def _set_kspace_coulomb(
    kspace_styles: list | None,
    force_calcs: ForceCalculations,
) -> None:
    """Set coulomb_type on force_calcs from kspace_style data."""
    if not kspace_styles:
        return
    kspace_entry = kspace_styles[0]
    if not (isinstance(kspace_entry, list) and len(kspace_entry) >= 1):
        return
    kspace_name = str(kspace_entry[0]).lower()
    coulomb_type = _KSPACE_COULOMB_TYPE.get(kspace_name)
    if coulomb_type is not None:
        force_calcs.coulomb_type = coulomb_type


def _set_neighbor_cutoff(
    neighbor: list | None,
    force_calcs: ForceCalculations,
    dist_unit: Any,
) -> None:
    """Set neighbor_update_cutoff on force_calcs from neighbor data."""
    if not neighbor:
        return
    nb_entry = neighbor[0]
    if not (isinstance(nb_entry, list) and len(nb_entry) >= 1):
        return
    try:
        skin = float(nb_entry[0]) * dist_unit
        existing = force_calcs.vdw_cutoff
        force_calcs.neighbor_update_cutoff = (
            existing + skin if existing is not None else skin
        )
    except (ValueError, TypeError):
        pass


def _set_neighbor_frequency(
    neigh_modify: list | None,
    force_calcs: ForceCalculations,
) -> None:
    """Set neighbor_update_frequency on force_calcs from neigh_modify data."""
    if not neigh_modify:
        return
    nm_entry = neigh_modify[0]
    if not isinstance(nm_entry, list):
        return
    try:
        every_idx = nm_entry.index('every')
        force_calcs.neighbor_update_frequency = int(nm_entry[every_idx + 1])
    except (ValueError, IndexError):
        pass


def _build_harmonic_bond(
    style: str,
    coeffs: list | None,
    distance_unit: Any,
    energy_unit: Any,
) -> HarmonicBond:
    """Create a HarmonicBond potential from LAMMPS harmonic bond coefficients."""
    potential = HarmonicBond()
    potential.type = 'bond'
    potential.functional_form = 'harmonic'
    potential.name = style
    if coeffs is not None and len(coeffs) > 0:
        try:
            n_interactions = len(coeffs)
            k_vals = np.array([float(c[1]) for c in coeffs])
            r0_vals = np.array([float(c[2]) for c in coeffs])
            potential.n_interactions = n_interactions
            potential.force_constant = k_vals.mean() * energy_unit / distance_unit**2
            potential.equilibrium_value = r0_vals.mean() * distance_unit
        except (IndexError, TypeError, ValueError):
            pass
    return potential


def _build_harmonic_angle(
    style: str,
    coeffs: list | None,
    energy_unit: Any,
) -> HarmonicAngle:
    """Create a HarmonicAngle potential from LAMMPS harmonic angle coefficients."""
    potential = HarmonicAngle()
    potential.type = 'angle'
    potential.functional_form = 'harmonic'
    potential.name = style
    if coeffs is not None and len(coeffs) > 0:
        try:
            n_interactions = len(coeffs)
            k_vals = np.array([float(c[1]) for c in coeffs])
            theta0_vals = np.array([float(c[2]) for c in coeffs])
            potential.n_interactions = n_interactions
            potential.force_constant = k_vals.mean() * energy_unit / _ureg.radian**2
            potential.equilibrium_value = (theta0_vals.mean() * _ureg.degree).to(
                'radian'
            )
        except (IndexError, TypeError, ValueError):
            pass
    return potential


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
