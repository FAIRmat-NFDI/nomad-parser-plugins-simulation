import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Program, Simulation
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

    def parse_system(self, simulation: Simulation) -> None:
        n_traj = self.traj_parsers.eval('n_frames')
        if n_traj is None:
            return
        self.n_atoms = [self.traj_parsers.eval('get_n_atoms', n) for n in range(n_traj)]
        self.trajectory_steps = [
            step
            for n in range(n_traj)
            if (step := self.traj_parsers.eval('get_step', n)) is not None
        ]

        for step in self.trajectory_steps:
            traj_n = self.trajectory_steps.index(step)
            lattice_vectors = self.traj_parsers.eval('get_lattice_vectors', traj_n)
            if lattice_vectors is not None:
                lattice_vectors = self.apply_unit(lattice_vectors, 'distance')
            velocities = self.traj_parsers.eval('get_velocities', traj_n)
            if velocities is not None:
                velocities = self.apply_unit(velocities, 'velocity')
            if traj_n == 0:  # TODO add references to the bond list for other steps
                # TODO: update get_bond_list_from_model_contributions,
                # TODO: maybe move to MDParserUtils?
                # bond_list = get_bond_list_from_model_contributions(
                #     sec_run, method_index=-1, model_index=-1
                # )
                if self._bond_list is None:
                    # Convert bond list returned by data parser from
                    # List[Tuple[None, np.ndarray]] to List[Tuple[int, int]]
                    bonds = self._data_parser.get('Bonds', None)
                    if bonds is None or bonds[0][1].size == 0:
                        self._bond_list = None
                    else:
                        self._bond_list = list(
                            map(tuple, bonds[0][1][:, 2:4].astype(int))
                        )
            # Set the structure of the data dictionary according to ModelSystem
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
            }
            self.parse_trajectory_step(particles_dict, simulation)

        # parse atomsgroup (moltypes --> molecules --> residues)
        # Only information from first frame is used to date
        first_frame = 0
        atoms_info = self._mdanalysistraj_parser.get('atoms_info', None)
        if atoms_info is None:
            atoms_info = self.traj_parsers.eval('atoms_info')
            if isinstance(atoms_info, list):
                atoms_info = (
                    atoms_info[first_frame] if atoms_info else None
                )  # using info from the initial frame
        if atoms_info is not None:
            # atoms_moltypes = np.array(atoms_info.get('moltypes', []))
            # atoms_molnums = np.array(atoms_info.get('molnums', []))
            # atoms_resids = np.array(atoms_info.get('resids', []))
            atoms_elements = np.array(
                atoms_info.get('elements', ['CGX'] * self.n_atoms)
            )
            atoms_types = np.array(atoms_info.get('types', []))
            atom_labels = [
                particle_state.label
                for particle_state in simulation.model_system[
                    first_frame
                ].particle_states
            ]
            if 'CGX' in atoms_elements:
                atoms_elements = (
                    np.array(atom_labels)
                    if atom_labels and 'CGX' not in atom_labels
                    else atoms_types
                )

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
            traj_file: str, index: int, file_type: str, data_file: str
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
