import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem
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

    def apply_unit(self, value: Any, unit: str) -> float:
        if not hasattr(value, 'units'):
            value = value * self._log_parser.units.get(unit, 1)
        return value

    def parse_method(self, simulation: Simulation) -> None:
        # TODO: replace with counterparts from nomad_simulations!

        if self.traj_parsers[0].mainfile is None or self._data_parser.mainfile is None:
            return

        if self.traj_parsers.eval('n_frames') is None:
            return

        masses = self._data_parser.get('Masses', None)
        self.traj_parsers[0].masses = masses

        # TODO: find best place for first attempt to set _bond_list
        # Extract bond list from MDAnalysis universe
        self._bond_list = [
            tuple(interaction['atom_indices'])
            for interaction in self._mdanalysistraj_parser.get_interactions()
            if interaction['type'] == 'bond'
        ]

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
                # Convert from List[Tuple[None, np.ndarray]] to List[Tuple[int, int]]
                self._bond_list = list(map(tuple, bonds[0][1][:, 2:4].astype(int)))

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
                'cell': {
                    'lattice_vectors': lattice_vectors,
                    'periodic_boundary_conditions': self.traj_parsers.eval(
                        'get_pbc', traj_n
                    ),
                },
                'labels': self.traj_parsers.eval('get_atom_labels', traj_n),
                'n_particles': self.traj_parsers.eval('get_n_atoms', traj_n),
                'positions': self.apply_unit(
                    self.traj_parsers.eval('get_positions', traj_n), 'distance'
                ),
                'velocities': velocities,
                'bond_list': self._bond_list if self._bond_list else None,
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
            is_representative=(i_molecule == 0),
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
            file_type = traj_file.split('.')[-1]

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
        # # include input controls from log file
        # self.parse_input()

        # # parse thermodynamic data from log file
        # self.parse_thermodynamic_data()

        # self.parse_workflow()

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
