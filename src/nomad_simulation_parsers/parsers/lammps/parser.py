import os
from collections.abc import Iterable
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_method import ModelMethod
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

        self.units = self._log_parser.units

    def apply_unit(self, value: Any, unit: str) -> float:
        if not hasattr(value, 'units'):
            value = value * self.units.get(unit, 1)
        return value

    def parse_method(self, simulation: Simulation) -> None:
        # TODO: replace with counterparts from nomad_simulations!

        if self.traj_parsers[0].mainfile is None or self._data_parser.mainfile is None:
            return

        if self.traj_parsers.eval('n_frames') is None:
            return

        sec_method = ModelMethod()
        # sec_force_field = ForceField()
        # sec_method.force_field = sec_force_field
        # sec_model = Model()
        # sec_force_field.model.append(sec_model)

        # Old parsing of method with text parser
        masses = self._data_parser.get('Masses', None)
        self.traj_parsers[0].masses = masses

        # ? Should _bond_list be set here, or is there a better place for it?
        # ? If anything failed with the universe, would it have errored out already, or happen here?
        # Exract bond list from MDAnalysis universe
        self._bond_list = [
            tuple(interaction['atom_indices'])
            for interaction in self._mdanalysistraj_parser.get_interactions()
            if interaction['type'] == 'bond'
        ]

        # parse method with MDAnalysis (should be a backup for the charges and masses...
        # but the interactions are most easily read from the MDA universe right now)
        # n_atoms = self.traj_parsers.eval('get_n_atoms', 0)
        # if n_atoms is not None:
        #     atoms_info = self._mdanalysistraj_parser.get('atoms_info', None)
        #     labels = self.traj_parsers.eval('labels')
        #     # TODO: LB, test
        #     # if labels is None or 'CGX' in labels:
        #     #     atom_types = self._mdanalysistraj_parser.get('types', None) # atom_types = self._mdanalysis.get('atom_types')
        #     #     #check if none and revert to X's if none
        #     #     if atom_types is None:
        #     #         atom_types = atoms_info.get('types', None)
        #     #         #atom_types = ['CGX']*n_atoms
        #     #     else:
        #     #         labels = [f'X_{atom_type}' for atom_type in atom_types]
        #     for n in range(n_atoms):
        #         sec_atom = AtomParameters()
        #         sec_method.atom_parameters.append(sec_atom)
        #         sec_atom.charge = atoms_info.get('charges', [None] * (n + 1))[n]
        #         sec_atom.mass = atoms_info.get('masses', [None] * (n + 1))[n]
        #         # TODO: LB edit, test
        #         # sec_atom.label = (labels[n] if labels is not None else f'X_{atom_types[n]}')  # *[n]

    #     # TODO address case types are numbered instead of giving atom labels (fix tests accordingly)
    #     interactions = self._mdanalysistraj_parser.get_interactions()
    #     # for interaction in interactions:
    #     #     for key, val in interaction.items():
    #     #         quantity_def = Interaction.m_def.all_quantities.get(key)
    #     #         if quantity_def and quantity_def.shape:
    #     #             # TODO reshape properly
    #     #             interaction[key] = [val]
    #     self.parse_interactions(interactions, sec_model)

    #     # Force Calculation Parameters
    #     sec_force_calculations = ForceCalculations()
    #     sec_force_field.force_calculations = sec_force_calculations
    #     for pairstyle in self._log_parser.get('pair_style', []):
    #         pairstyle_args = pairstyle[1:]
    #         pairstyle = pairstyle[0].lower()
    #         if (
    #             'lj' in pairstyle and 'coul' not in pairstyle
    #         ):  # only cover the simplest case
    #             sec_force_calculations.vdw_cutoff = (
    #                 float(pairstyle_args[-1]) * ureg.nanometer
    #                 # TODO: LB edit
    #                 # float(pairstyle_args[-1]) * ureg.angstrom
    #             )
    #         if 'coul' in pairstyle:
    #             if 'streitz' in pairstyle:
    #                 cutoff = float(pairstyle_args[0])
    #             else:
    #                 cutoff = float(pairstyle_args[-1])
    #             sec_force_calculations.coulomb_cutoff = cutoff * ureg.nanometer
    #             # TODO: LB edit
    #             # sec_force_calculations.coulomb_cutoff = cutoff * ureg.angstrom
    #         val = self._log_parser.get('kspace_style', None)
    #         if val is not None:
    #             kspacestyle = val[0][0].lower()
    #             if 'ewald' in kspacestyle:
    #                 sec_force_calculations.coulomb_type = 'ewald'
    #             elif 'pppm' in kspacestyle:
    #                 sec_force_calculations.coulomb_type = (
    #                     'particle_particle_particle_mesh'
    #                 )
    #             elif 'msm' in kspacestyle:
    #                 sec_force_calculations.coulomb_type = 'multilevel_summation'

    #     sec_neighbor_searching = NeighborSearching()
    #     sec_force_calculations.neighbor_searching = sec_neighbor_searching
    #     val = self._log_parser.get('neighbor', None)
    #     if val is not None:
    #         neighbor = val[0][0]  # just use the first instance for now
    #         vdw_cutoff = sec_force_calculations.vdw_cutoff
    #         if vdw_cutoff is not None:
    #             sec_neighbor_searching.neighbor_update_cutoff = (
    #                 float(neighbor) * ureg.nanometer
    #             )
    #             sec_neighbor_searching.neighbor_update_cutoff += vdw_cutoff
    #     val = self._log_parser.get('neigh_modify', None)
    #     if val is not None:
    #         neighmodify = val[0]  # just use the first instace for now
    #         neighmodify = np.array([str(i).lower() for i in neighmodify])
    #         if 'every' in neighmodify:
    #             index = np.where(neighmodify == 'every')[0]
    #             sec_neighbor_searching.neighbor_update_frequency = int(
    #                 neighmodify[index + 1]
    #            )

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
                # TODO: update get_bond_list_from_model_contributions, maybe move to MDParserUtils?
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
        # ! Only information from first frame is used to date
        first_frame = 0
        atoms_info = self._mdanalysistraj_parser.get('atoms_info', None)
        if atoms_info is None:
            atoms_info = self.traj_parsers.eval('atoms_info')
            if isinstance(atoms_info, list):
                atoms_info = (
                    atoms_info[first_frame] if atoms_info else None
                )  # using info from the initial frame
        if atoms_info is not None:
            atoms_moltypes = np.array(atoms_info.get('moltypes', []))
            atoms_molnums = np.array(atoms_info.get('molnums', []))
            atoms_resids = np.array(atoms_info.get('resids', []))
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
            atoms_resnames = np.array(atoms_info.get('resnames', []))
            moltypes = np.unique(atoms_moltypes)
            # for i_moltype, moltype in enumerate(moltypes):
            #     # Only add atomsgroup for initial system for now
            #     # ! AtomsGroup deprecated, sub_system = ModelSystem() now!
            #     sec_molecule_group = ModelSystem()
            #     # TODO: append to simulation?
            #     # sec_run.system[0].atoms_group.append(sec_molecule_group)
            #     # ? ModelSystem has no 'label', 'ParticleState' does.
            #     # ? Using a 'ParticleState' class for a group of molecules is counter-intuitive
            #     sec_molecule_group.branch_label = f'group_{moltype}'
            #     # ? 'molecule_group' deprecated,
            #     # ? which of 'cluster', 'bulk', or 'unavailable' to use?
            #     sec_molecule_group.type = 'cluster'  # 'molecule_group'
            #     sec_molecule_group.index = i_moltype
            #     sec_molecule_group.atom_indices = np.where(atoms_moltypes == moltype)[0]
            #     sec_molecule_group.n_atoms = len(sec_molecule_group.atom_indices)
            #     sec_molecule_group.is_molecule = False
            #     # mol_nums is the molecule identifier for each atom
            #     mol_nums = atoms_molnums[sec_molecule_group.atom_indices]
            #     moltype_count = np.unique(mol_nums).shape[0]
            #     sec_molecule_group.composition_formula = f'{moltype}({moltype_count})'

            #     molecules = atoms_molnums
            #     for i_molecule, molecule in enumerate(
            #         np.unique(molecules[sec_molecule_group.atom_indices])
            #     ):
            #         sec_molecule = ModelSystem()
            #         # TODO: append to sec_molecule_group
            #         # sec_molecule_group.atoms_group.append(sec_molecule)
            #         if i_molecule == 0:
            #             sec_molecule.is_representative = True
            #         sec_molecule.index = i_molecule
            #         sec_molecule.atom_indices = np.where(molecules == molecule)[0]
            #         sec_molecule.n_atoms = len(sec_molecule.atom_indices)
            #         # use first particle to get the moltype
            #         # not sure why but this value is being cast to int, cast back to str
            #         sec_molecule.label = str(
            #             atoms_moltypes[sec_molecule.atom_indices[0]]
            #         )
            #         sec_molecule.type = 'molecule'
            #         sec_molecule.is_molecule = True

            #         mol_resids = np.unique(atoms_resids[sec_molecule.atom_indices])
            #         n_res = mol_resids.shape[0]
            #         if n_res == 1:
            #             elements = atoms_elements[sec_molecule.atom_indices]
            #             sec_molecule.composition_formula = get_composition(elements)
            #         else:
            #             mol_resnames = atoms_resnames[sec_molecule.atom_indices]
            #             restypes = np.unique(mol_resnames)
            #             for i_restype, restype in enumerate(restypes):
            #                 sec_monomer_group = ModelSystem()
            #                 # TODO: append to sec_molecule
            #                 # sec_molecule.atoms_group.append(sec_monomer_group)
            #                 restype_indices = np.where(atoms_resnames == restype)[0]
            #                 sec_monomer_group.label = f'group_{restype}'
            #                 sec_monomer_group.type = 'monomer_group'
            #                 sec_monomer_group.index = i_restype
            #                 sec_monomer_group.atom_indices = np.intersect1d(
            #                     restype_indices, sec_molecule.atom_indices
            #                 )
            #                 sec_monomer_group.n_atoms = len(
            #                     sec_monomer_group.atom_indices
            #                 )
            #                 sec_monomer_group.is_molecule = False

            #                 restype_resids = np.unique(
            #                     atoms_resids[sec_monomer_group.atom_indices]
            #                 )
            #                 restype_count = restype_resids.shape[0]
            #                 sec_monomer_group.composition_formula = (
            #                     f'{restype}({restype_count})'
            #                 )
            #                 for i_res, res_id in enumerate(restype_resids):
            #                     sec_residue = ModelSystem()
            #                     # TODO: append to sec_monomer_group
            #                     # sec_monomer_group.atoms_group.append(sec_residue)
            #                     sec_residue.index = i_res
            #                     atom_indices = np.where(atoms_resids == res_id)[0]
            #                     sec_residue.atom_indices = np.intersect1d(
            #                         atom_indices, sec_monomer_group.atom_indices
            #                     )
            #                     sec_residue.n_atoms = len(sec_residue.atom_indices)
            #                     sec_residue.label = str(restype)
            #                     sec_residue.type = 'monomer'
            #                     sec_residue.is_molecule = False
            #                     elements = atoms_elements[sec_residue.atom_indices]
            #                     sec_residue.composition_formula = get_composition(
            #                         elements
            #                     )

            #             names = atoms_resnames[sec_molecule.atom_indices]
            #             ids = atoms_resids[sec_molecule.atom_indices]
            #             # filter for the first instance of each residue, as to not overcount
            #             __, ids_count = np.unique(ids, return_counts=True)
            #             # get the index of the first atom of each residue
            #             ids_firstatom = np.cumsum(ids_count)[:-1]
            #             # add the 0th index manually
            #             ids_firstatom = np.insert(ids_firstatom, 0, 0)
            #             names_firstatom = names[ids_firstatom]
            #             sec_molecule.composition_formula = get_composition(
            #                 names_firstatom
            #             )

    def write_to_archive(self) -> None:
        self.archive.data = Simulation(program=Program(name='LAMMPS'))
        # LAMMPS mainfile is the main log file
        self.basename = os.path.basename(self.mainfile)
        self.basedir = os.path.dirname(self.mainfile)
        self._log_parser.mainfile = self.mainfile
        self._log_parser.logger = self.logger
        self._log_parser._units = None

        # parse data from auxiliary log file
        if self._log_parser.get('log') is not None:
            self._aux_log_parser.mainfile = os.path.join(
                self._log_parser.maindir,
                self._log_parser.get('log')[0],
            )
            # we assign units here which is read from log parser
            self._aux_log_parser._units = self._log_parser.units
            self._aux_log_parser.logger = self.logger

        self._traj_parser.logger = self.logger
        self._traj_parser._chemical_symbols = None
        self._xyztraj_parser.logger = self.logger
        self._mdanalysistraj_parser.logger = self.logger
        # self._mdparser = MDParser()
        # self._mdparser.logger = self.logger
        self._data_parser.logger = self.logger

        # parse data file associated with calculation
        data_files = self._log_parser.get_data_files()
        if len(data_files) > 1:
            self.logger.warning('Multiple data files are specified')
        if data_files:
            self._data_parser.mainfile = data_files[0]

        # parse trajectorty file associated with calculation
        traj_files = self._log_parser.get_traj_files()
        if len(traj_files) > 1:
            self.logger.warning('Multiple traj files are specified')

        parsers = []
        for n, traj_file in enumerate(traj_files):
            # parser initialization for each traj file cannot be avoided as there are
            # cases where traj files can share the same parser
            file_type = self._log_parser.get(
                'dump', [[1, 'all', traj_file.split('.')[-1]]] * (n + 1)
            )[n][2]
            # TODO: add support for other LAMMPs dump file formats (https://docs.lammps.org/dump.html)
            if file_type == 'dcd' and data_files:
                traj_parser = MDAnalysisParser(topology_format='DATA', format='DCD')
                traj_parser.mainfile = data_files[0]
                traj_parser.auxilliary_files = [traj_file]
                self._mdanalysistraj_parser = traj_parser
            elif file_type == 'xyz' and data_files:
                traj_parser = MDAnalysisParser(topology_format='DATA', format='XYZ')
                traj_parser.mainfile = data_files[0]
                traj_parser.auxilliary_files = [traj_file]
                self._mdanalysistraj_parser = traj_parser
            # TODO: LB edit
            # elif file_type == 'atom' and data_files:
            #     traj_parser = MDAnalysisParser(
            #         topology_format='DATA', format='LAMMPSDUMP'
            #     )
            #     if data_files:
            #         traj_parser.mainfile = data_files[0]
            #     traj_parser.auxilliary_files = [traj_file]

            #     if traj_parser.universe is None or 'CGX' in traj_parser.get(
            #         'atoms_info', {}
            #     ).get('names', []):
            #         # mda necessary to calculate rdf and atomsgroup
            #         if n == 0:
            #             self._mdanalysistraj_parser = traj_parser
            #         traj_parser = TrajParser()
            #         traj_parser.mainfile = traj_file
            elif file_type == 'custom' and data_files:
                custom_options = self._log_parser.get('dump')[n][5:]
                custom_options = [
                    option.replace('xu', 'x') for option in custom_options
                ]
                custom_options = [
                    option.replace('yu', 'y') for option in custom_options
                ]
                custom_options = [
                    option.replace('zu', 'z') for option in custom_options
                ]
                custom_options = ' '.join(custom_options)
                traj_parser = MDAnalysisParser(
                    topology_format='DATA',
                    format='LAMMPSDUMP',
                    atom_style=custom_options,
                )
                if data_files:
                    traj_parser.mainfile = data_files[0]
                traj_parser.auxilliary_files = [traj_file]
                # try to check if MDAnalysis can construct the universe or at least parse
                # the atoms, otherwise will fall back to TrajParser
                if traj_parser.universe is None or 'CGX' in traj_parser.get(
                    'atoms_info', {}
                ).get('names', []):
                    # mda necessary to calculate rdf and atomsgroup
                    if n == 0:
                        self._mdanalysistraj_parser = traj_parser
                    traj_parser = TrajParser()
                    traj_parser.mainfile = traj_file
            else:
                self.logger.warning(f'File type of {traj_file} not recognized.')
                traj_parser = TrajParser()
                traj_parser.mainfile = traj_file
                # TODO provide support for other file types
            parsers.append(traj_parser)

        self.traj_parsers = TrajParsers(parsers)
        if self.traj_parsers[0] is None:
            return
        # parse data from auxiliary log file
        if self._log_parser.get('log') is not None:
            self._aux_log_parser.mainfile = os.path.join(
                self._log_parser.maindir, self._log_parser.get('log')[0]
            )
            # we assign units here which is read from log parser
            self._aux_log_parser._units = self._log_parser.units

        self.parse_method(self.archive.data)

        self.parse_system(self.archive.data)

        # # include input controls from log file
        # self.parse_input()

        # # parse thermodynamic data from log file
        # self.parse_thermodynamic_data()

        # self.parse_workflow()

        self._mdanalysistraj_parser.close()
        for parser in parsers:
            parser.close()


class LammpsParser(MatchingParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.archive_writer = LammpsArchiveWriter()

    # ? Really needed for the LAMMPS parser?
    # ? Would it make sense to handle a potential auxillary log file here,
    # ? since there seems to be an issue with its handling in the old parser?
    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable[str]:
        """
        TODO: Documentation
        """
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )

        if is_mainfile:
            # ? Handle check for auxiliary log file here?
            return is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
