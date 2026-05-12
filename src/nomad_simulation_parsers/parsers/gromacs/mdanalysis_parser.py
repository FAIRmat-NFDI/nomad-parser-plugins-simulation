import operator
import re
from collections.abc import Callable
from typing import Any

import MDAnalysis
import numpy as np
from ase import Atoms
from ase.data import chemical_symbols
from MDAnalysis.topology.tpr import setting as tpr_setting
from MDAnalysis.topology.tpr import utils as tpr_utils
from nomad.normalizing.optimade import normalized_atom_labels

from nomad_simulation_parsers.parsers.utils.mdanalysisparser import MDAnalysisParser

# =============================================================================
# GROMACS TPR PARSING WITH MDANALYSIS
# =============================================================================
#
# WHAT MDANALYSIS PROVIDES:
# -------------------------
# 1. Topology connectivity:
#    - universe.bonds: Bond pairs (atom indices)
#    - universe.angles: Angle triplets
#    - universe.dihedrals: Dihedral quadruplets
#    - universe.impropers: Improper dihedral quadruplets
#
# 2. Atom properties:
#    - Atom names, types, residues
#    - Charges, masses
#
# 3. Coordinates:
#    - Positions, velocities, forces (if present in TPR)
#
# 4. System properties:
#    - Box dimensions
#    - Number of atoms
#
# WHAT MDANALYSIS DOES NOT PROVIDE:
# ---------------------------------
# 1. Force field parameters:
#    - Bond force constants and equilibrium lengths
#    - Angle force constants and equilibrium angles
#    - Dihedral parameters
#    - Lennard-Jones σ, ε values
#
# 2. Parameter-to-topology mapping:
#    - Which bond uses which parameter set
#    - Parameter set indices from ilist section
#
# 3. Numerical settings:
#    - Cutoff distances (read from mdp/log instead)
#    - Neighbor list parameters
#    - PME grid settings
#
# CUSTOM IMPLEMENTATION BELOW:
# ---------------------------
# - get_force_field_parameters(): Reads parameter VALUES from TPR
#   Uses MDAnalysis TPXUnpacker to access binary data directly
#   Returns parameter arrays but cannot connect them to specific interactions
#
# - get_interactions(): Inherited from MDAnalysisParser
#   Provides topology (atom indices) without parameters
#
# For full force field support, would need to:
# 1. Parse ilist section (not exposed by MDAnalysis)
# 2. Match each interaction to its parameter set index
# 3. Combine topology + parameters into complete ForceField objects
# =============================================================================


class GromacsMDAnalysisParser(MDAnalysisParser):
    def reset(self):
        super().reset()
        self.data = None
        self._func_unpacker = {}
        self.data = None

    @property
    def func_unpacker(self) -> dict[int, Callable]:
        if not self.data:
            return self._func_unpacker
        if not self._func_unpacker:
            self._func_unpacker = {
                tpr_setting.F_ANGLES: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_G96ANGLES: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_BONDS: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_G96BONDS: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_HARMONIC: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_IDIHS: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_RESTRANGLES: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_LINEAR_ANGLES: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_FENEBONDS: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_RESTRBONDS: self.eval_unpacker(['R'] * 8),
                tpr_setting.F_TABBONDS: self.eval_unpacker(['R', 'I', 'R']),
                tpr_setting.F_TABBONDSNC: self.eval_unpacker(['R', 'I', 'R']),
                tpr_setting.F_TABANGLES: self.eval_unpacker(['R', 'I', 'R']),
                tpr_setting.F_TABDIHS: self.eval_unpacker(['R', 'I', 'R']),
                tpr_setting.F_CROSS_BOND_BONDS: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_CROSS_BOND_ANGLES: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_UREY_BRADLEY: self.eval_unpacker(
                    ['R'] * (8 if self.check_header(79) else 4)
                ),
                tpr_setting.F_QUARTIC_ANGLES: self.eval_unpacker(['R', '5R']),
                tpr_setting.F_BHAM: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_MORSE: self.eval_unpacker(
                    ['R'] * (6 if self.check_header(79) else 3)
                ),
                tpr_setting.F_CUBICBONDS: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_POLARIZATION: self.eval_unpacker(['R']),
                tpr_setting.F_ANHARM_POL: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_WATER_POL: self.eval_unpacker(['R'] * 6),
                tpr_setting.F_THOLE_POL: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_LJ: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_LJ14: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_LJC14_Q: self.eval_unpacker(['R'] * 5),
                tpr_setting.F_LJC_PAIRS_NB: self.eval_unpacker(['R'] * 4),
                tpr_setting.F_PIDIHS: self.eval_unpacker(['R'] * 4 + ['I']),
                tpr_setting.F_ANGRES: self.eval_unpacker(['R'] * 4 + ['I']),
                tpr_setting.F_ANGRESZ: self.eval_unpacker(['R'] * 4 + ['I']),
                tpr_setting.F_PDIHS: self.eval_unpacker(['R'] * 4 + ['I']),
                tpr_setting.F_RESTRDIHS: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_DISRES: self.eval_unpacker(['I'] * 2 + ['R'] * 4),
                tpr_setting.F_ORIRES: self.eval_unpacker(['I'] * 3 + ['R'] * 3),
                tpr_setting.F_DIHRES: self.eval_unpacker(
                    ['I'] * 2 + ['R'] * (6 if self.check_header(72) else 3)
                ),
                tpr_setting.F_POSRES: self.eval_unpacker(['RVEC'] * 4),
                tpr_setting.F_FBPOSRES: self.eval_unpacker(['I', 'RVEC', 'R', 'R']),
                tpr_setting.F_CBTDIHS: self.eval_unpacker(
                    [f'{tpr_setting.NR_CBTDIHS}R']
                ),
                tpr_setting.F_RBDIHS: self.eval_unpacker(
                    [f'{tpr_setting.NR_RBDIHS}R'] * 2
                ),
                tpr_setting.F_FOURDIHS: self.eval_unpacker(
                    [f'{tpr_setting.NR_RBDIHS}R'] * 2
                ),
                tpr_setting.F_CONSTR: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_CONSTRNC: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_SETTLE: self.eval_unpacker(['R'] * 2),
                tpr_setting.F_VSITE1: self.eval_unpacker([]),
                tpr_setting.F_VSITE2: self.eval_unpacker(['R']),
                tpr_setting.F_VSITE2FD: self.eval_unpacker(['R']),
                tpr_setting.F_VSITE3: self.eval_unpacker(['RFL']),
                tpr_setting.F_VSITE3FD: self.eval_unpacker(['RFL']),
                tpr_setting.F_VSITE3FAD: self.eval_unpacker(['RFL']),
                tpr_setting.F_VSITE3OUT: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_VSITE4FD: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_VSITE4FDN: self.eval_unpacker(['R'] * 3),
                tpr_setting.F_VSITEN: self.eval_unpacker(['I', 'R']),
                tpr_setting.F_GB12: self.eval_unpacker(
                    ['R'] * (9 if self.check_header(68, op=operator.lt) else 5)
                ),
                tpr_setting.F_GB13: self.eval_unpacker(
                    ['R'] * (9 if self.check_header(68, op=operator.lt) else 5)
                ),
                tpr_setting.F_GB14: self.eval_unpacker(
                    ['R'] * (9 if self.check_header(68, op=operator.lt) else 5)
                ),
                tpr_setting.F_CMAP: self.eval_unpacker(['I'] * 2),
            }

        return self._func_unpacker

    def check_header(self, value: int, key='fver', op=operator.ge) -> bool | None:
        if self.data is None:
            return None
        return op(getattr(self.header, key), value)

    def eval_unpacker(self, unpackers: list[str]) -> Callable:
        def func() -> list[Any]:
            parameters = []
            for name in unpackers:
                if name == 'I':
                    p = self.data.unpack_int()
                elif name == 'R':
                    p = self.data.unpack_real()
                elif name == 'RVEC':
                    p = tpr_utils.ndo_rvec(self.data)
                elif name == 'RFL':
                    p = self.data.unpack_reafilel()
                elif m := re.match(r'(\d+)R', name):
                    p = tpr_utils.ndo_real(self.data, int(m.group(1)))
                else:
                    self.logger.error('Unrecognized unpack method.')
                    p = None
                parameters.append(p)
            return parameters

        return func

    def get_interactions(self, gromacs_version: str = None) -> list[dict[str, Any]]:
        interactions = super().get_interactions()
        # add force field parameters
        try:
            interactions.extend(self.get_force_field_parameters(gromacs_version))
        except Exception:
            self.logger.warning('Error parsing force field parameters.')

        self._results['interactions'] = interactions

        return interactions

    def get_force_field_parameters(
        self, gromacs_version: str = None
    ) -> list[dict[str, Any]]:
        """
        Read force field parameters from GROMACS TPR binary file.

        TPR FILE STRUCTURE:
        ==================
        The TPR file stores force field parameters in a compact binary format:

        1. HEADER: Version info, system size, flags (bTop, bBox, bX, bV, bF)
        2. TOPOLOGY SECTION (if bTop == True):
           - Symbol table (atom/residue names)
           - ntypes: number of parameter sets
           - functypes[ntypes]: array of function type IDs
           - reppow: repulsion power (double)
           - fudgeQQ: 1-4 Coulomb scaling factor (real)
           - For each functype:
             * Read N parameters based on type specification
             * Store as {'type': name, 'parameters': [...]}

        PARAMETER TYPES:
        ===============
        Each functype ID determines the interaction type and parameter count:

        F_BONDS (0):          4 reals - Harmonic bond (k, r0, ?, ?)
        F_G96BONDS (1):       4 reals - GROMOS96 bond
        F_MORSE (2):          2-6 reals - Morse potential (D, alpha, r0, ...)
        F_ANGLES (10):        4 reals - Harmonic angle (k, θ0, ?, ?)
        F_PDIHS (19):         4R + 1I - Proper dihedral (φs, kφ, mult)
        F_RBDIHS (27):        12R - Ryckaert-Bellemans dihedral
        F_LJ (37):            2 reals - Lennard-Jones (C6, C12 or σ, ε)
        F_LJ14 (45):          4 reals - LJ 1-4 interactions
        F_CONSTR (62):        2 reals - Bond constraint (b0, tolerance)
        F_SETTLE (64):        2 reals - Water constraint (dOH, dHH)

        CURRENT LIMITATION:
        ==================
        This implementation extracts parameter VALUES only, not their assignments.
        The TPR file contains:
        - Parameter sets (this function extracts these)
        - Topology lists (bonds, angles, dihedrals) with atom indices + functype index

        For full force field parsing, we would need to:
        1. Parse interaction lists (ilist section after topology)
        2. Match each interaction (e.g., bond 0-1) to its parameter set
        3. Create proper ForceField.contributions with:
           - functional_form (e.g., 'harmonic_bond')
           - particle_indices from topology
           - parameters from this function
        4. Handle version-dependent formats and unit conversions

        The interaction lists are stored later in the TPR file structure:
        - ilist[F_BONDS]: list of [atom_i, atom_j, param_index]
        - ilist[F_ANGLES]: list of [atom_i, atom_j, atom_k, param_index]
        - etc.

        Currently, MDAnalysis reads the topology and provides interactions via
        universe.bonds, universe.angles, etc., but does NOT connect them to the
        parameter values extracted here. That would require extending MDAnalysis
        or implementing a custom TPR reader.

        RETURNS:
        =======
        List of dicts: [{'type': 'LJ (SR)', 'parameters': [C6, C12]}, ...]
        """
        # copied from MDAnalysis.topology.tpr.utils
        # TODO Revamp interactions section to only extract meaningful info
        if MDAnalysis.__version__.split('.')[0] != '2':
            self.logger.warning(
                'MDAnalysis >= 2.0.0 is required for reading force field from tpr.'
                'Interactions will not be stored'
            )
            return []
        if not self.universe:
            if isinstance(self.universe_error, NotImplementedError):
                self.logger.warning(
                    'GROMACS TPR file version currently not supported by MDAnalysis. '
                    'Force field interactions will not be stored.'
                )
            elif self.universe_error:
                self.logger.warning(
                    'Failed to read TPR file: %s. '
                    'Force field interactions will not be stored.',
                    str(self.universe_error),
                )
            return []

        with open(self.mainfile, 'rb') as f:
            self.data = tpr_utils.TPXUnpacker(f.read())

        interactions: list[dict[str, Any]] = []

        # read header
        self.header = tpr_utils.read_tpxheader(self.data)
        # address compatibility issue
        if self.header.fver >= tpr_setting.tpxv_AddSizeField and self.check_header(
            27, 'fgen'
        ):
            actual_body_size = len(self.data.get_buffer()) - self.data.get_position()
            if actual_body_size == 4 * self.header.sizeOfTprBody:
                self.logger.error('Unsupported tpr format.')
                return interactions
            self.data = tpr_utils.TPXUnpacker2020.from_unpacker(self.data)

        # read other unimportant parts
        if self.header.bBox:
            tpr_utils.extract_box_info(self.data, self.header.fver)
        if self.header.ngtc > 0:
            if self.check_header(60, op=operator.lt):
                tpr_utils.ndo_real(self.data, self.header.ngtc)
            tpr_utils.ndo_real(self.data, self.header.ngtc)
        if not self.header.bTop:
            return interactions

        tpr_utils.do_symstr(self.data, tpr_utils.do_symtab(self.data))
        self.data.unpack_int()
        ntypes = self.data.unpack_int()
        # functional types
        functypes = tpr_utils.ndo_int(self.data, ntypes)
        self.data.unpack_double() if self.check_header(66) else 12.0
        self.data.unpack_real()
        # read the ffparams
        for i in functypes:
            unpacker = self.func_unpacker.get(i)
            if unpacker is None:
                self.logger.error('Unknown force field functype.')
                continue
            interactions.append(
                dict(type=tpr_setting.interaction_types[i][1], parameters=unpacker())
            )

        return interactions

    def _parse_element_from_name(self, name: str) -> str:
        """Extract element symbol from atom name like 'H1', 'CA', etc."""
        # Try to match 1 or 2 letter element at start of name
        match = re.match(r'([A-Z][a-z]?)', name)
        if match:
            symbol = match.group(1)
            # Validate it's a real element
            if symbol in chemical_symbols[1:]:
                return symbol
        return 'CGX'  # Fallback for coarse-grained or unknown

    def _generate_chemical_formula(
        self, particle_indices: np.ndarray, particle_arrays: dict
    ) -> str:
        """Generate chemical formula from particle indices.

        Uses ASE for atomistic systems. Falls back to counting residue names
        for coarse-grained particles with unrecognized element symbols.
        """
        elements = particle_arrays['elements'][particle_indices]
        valid_symbols = set(chemical_symbols[1:])
        valid_mask = np.array([str(e) in valid_symbols for e in elements])

        if valid_mask.all():
            atoms = Atoms(symbols=list(elements))
            result = atoms.get_chemical_formula(mode='hill')
        else:
            resnames = particle_arrays['resnames'][particle_indices]
            unique_names, counts = np.unique(resnames, return_counts=True)
            result = ''.join(
                f'{name}{count}' if count > 1 else str(name)
                for name, count in zip(unique_names, counts)
            )
        return result

    def _get_particles_info_from_universe(self) -> dict[str, Any]:
        _atoms = self.universe.atoms
        particles_elements = np.array(
            [self._parse_element_from_name(atom.name) for atom in _atoms]
        )
        particles_types = np.array([atom.type for atom in _atoms])
        n_atoms = len(_atoms)
        moltypes = np.empty(n_atoms, dtype=object)
        molnums = np.zeros(n_atoms, dtype=int)
        segments = sorted(set([atom.segid for atom in _atoms]))
        for seg_idx, segment in enumerate(segments):
            atom_indices = np.array(
                [atom.index for atom in _atoms if atom.segid == segment]
            )
            n_residues = len(
                np.unique([atom.resid for atom in _atoms if atom.segid == segment])
            )

            # Determine molecule type identifier
            if n_residues == 1:
                # Single-residue molecule: use residue name
                moltype = _atoms[atom_indices[0]].resname
            else:
                # Multi-residue molecule: use segment ID
                moltype = segment

            moltypes[atom_indices] = moltype
            molnums[atom_indices] = seg_idx

        return {
            'moltypes': moltypes,
            'molnums': molnums,
            'resids': np.array([atom.resid for atom in _atoms]),
            'resnames': np.array([atom.resname for atom in _atoms]),
            'labels': [atom.name for atom in _atoms],
            'elements': particles_elements,
            'types': particles_types,
        }

    def _create_system_node(
        self, name: str | int, branch_label: str, particle_indices: np.ndarray, **kwargs
    ) -> dict[str, Any]:
        """
        Create a ModelSystem node with common setup.

        Args:
            name: System name
            branch_label: Hierarchy level label
            particle_indices: Indices of particles in this system
            **kwargs: Additional attributes: composition_formula, is_representative, ...
        """
        system = dict()
        system['name'] = str(name)
        system['branch_label'] = branch_label
        system['particle_indices'] = particle_indices
        system['is_molecule'] = branch_label == 'molecule'

        # Set any additional attributes
        for key, value in kwargs.items():
            system[key] = value

        return system

    def _create_molecule(
        self, molecule: int, i_molecule: int, particle_arrays: dict
    ) -> dict[str, Any]:
        """Create a single molecule with its residues."""

        def _create_residue(
            res_id: int,
            restype: str,
            parent_system: dict[str, Any],
            particle_arrays: dict,
        ) -> dict[str, Any]:
            """Create a single residue."""
            particle_indices = np.where(particle_arrays['resids'] == res_id)[0]
            particle_indices = np.intersect1d(
                particle_indices, parent_system['particle_indices']
            )

            composition_formula = self._generate_chemical_formula(
                particle_indices, particle_arrays
            )

            return self._create_system_node(
                name=f'{restype}',
                branch_label='monomer',
                particle_indices=particle_indices,
                composition_formula=composition_formula,
            )

        def _create_monomer_group(
            restype: str, parent_system: dict[str, Any], particle_arrays: dict[str, Any]
        ) -> dict[str, Any]:
            """Create a monomer group with its constituent residues."""

            restype_indices = np.where(particle_arrays['resnames'] == restype)[0]
            particle_indices = np.intersect1d(
                restype_indices, parent_system['particle_indices']
            )

            monomer_group = self._create_system_node(
                name=f'group_{restype}',
                branch_label='monomer_group',
                particle_indices=particle_indices,
            )

            # Add individual residues
            restype_resids = np.unique(
                particle_arrays['resids'][monomer_group['particle_indices']]
            )
            for res_id in restype_resids:
                residue = _create_residue(
                    res_id, restype, monomer_group, particle_arrays
                )
                monomer_group.setdefault('sub_systems', []).append(residue)

            # Set composition formula as RESTYPE(count)
            monomer_count = len(restype_resids)
            monomer_group['composition_formula'] = f'{restype}({monomer_count})'

            return monomer_group

        def _add_residue_hierarchy(
            sec_molecule: dict[str, Any], particle_arrays: dict[str, Any]
        ) -> None:
            """Add residue/monomer hierarchy to a molecule."""
            mol_resnames = particle_arrays['resnames'][sec_molecule['particle_indices']]
            restypes = np.unique(mol_resnames)

            for restype in restypes:
                sec_monomer_group = _create_monomer_group(
                    restype, sec_molecule, particle_arrays
                )
                sec_molecule.setdefault('sub_systems', []).append(sec_monomer_group)

        particle_indices = np.where(particle_arrays['molnums'] == molecule)[0]

        # Get molecule type from first atom in molecule
        moltype = particle_arrays['moltypes'][particle_indices[0]]
        molecule_name = f'{moltype}'

        mol_system = self._create_system_node(
            name=molecule_name,
            branch_label='molecule',
            particle_indices=particle_indices,
        )

        # Check if molecule has multiple residues
        mol_resids = np.unique(
            particle_arrays['resids'][mol_system['particle_indices']]
        )
        if len(mol_resids) > 1:
            _add_residue_hierarchy(mol_system, particle_arrays)
        else:
            # Single-residue molecule: add composition formula (terminal branch)
            mol_system['composition_formula'] = self._generate_chemical_formula(
                mol_system['particle_indices'], particle_arrays
            )

        return mol_system

    def parse_molecular_hierarchy(self) -> dict[str, Any]:
        """
        Build hierarchical subsystem structure from MDAnalysis universe.

        Creates up to 4 levels based on available topology:
        - Level 1: molecule_group (molecules of same type grouped)
        - Level 2: molecule (individual molecular fragment)
        - Level 3: monomer_group (residues of same type, for polymers)
        - Level 4: monomer (individual residue)

        Uses MDAnalysis residue information to build hierarchy. Computes composition
        formulas by counting atom types or residue types at each level. Returnes
        chemical formula for terminal branches (monomers or single-residue molecules).

        Args:
            source: Source dictionary (unused - we access MDAnalysis directly)
            **kwargs: Additional keyword arguments from annotation

        Returns:
            dict: Level-1 subsystems with nested sub_systems
        """

        mol_hierarchy: dict[str, Any] = {}

        def _get_particles_info() -> dict | None:
            """Get particle information from the first frame."""

            particles_info = self.get('atoms_info', None)

            if particles_info is None:
                particles_info = self._get_particles_info_from_universe()

            return particles_info

        def _extract_particle_arrays(particles_info: dict) -> dict:
            """Extract and process particle information arrays."""
            particle_labels = particles_info.get('labels', [])
            n_atoms = len(particle_labels)
            particles_elements = np.array(
                particles_info.get('elements', ['CGX'] * n_atoms), dtype=object
            )
            particles_types = np.array(particles_info.get('types', []))

            # Replace CGX placeholder elements per-atom using normalized_atom_labels,
            # which matches the longest valid chemical symbol found anywhere in the
            # atom name (e.g. OW→O, HW1→H, CA→Ca). Raw labels/types must not be
            # used directly since atom names and force-field types are not valid
            # ASE chemical symbols.
            cgx_mask = particles_elements == 'CGX'
            if cgx_mask.any():
                cgx_indices = np.where(cgx_mask)[0]
                fallback_names = [
                    particle_labels[i] if i < len(particle_labels) else ''
                    for i in cgx_indices
                ]
                particles_elements[cgx_indices] = normalized_atom_labels(fallback_names)

            return {
                'moltypes': np.array(particles_info.get('moltypes', [])),
                'molnums': np.array(particles_info.get('molnums', [])),
                'resids': np.array(particles_info.get('resids', [])),
                'resnames': np.array(particles_info.get('resnames', [])),
                'elements': particles_elements,
                'types': particles_types,
            }

        def _create_molecule_group(
            moltype: str, particle_arrays: dict
        ) -> dict[str, Any]:
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
                np.unique(molecules[molecule_group['particle_indices']])
            ):
                mol = self._create_molecule(molecule, i_molecule, particle_arrays)
                molecule_group.setdefault('sub_systems', []).append(mol)

            return molecule_group

        particles_info = _get_particles_info()
        if particles_info is None:
            return {}
        particle_arrays = _extract_particle_arrays(particles_info)

        # Build molecular hierarchy
        moltypes = np.unique(particle_arrays['moltypes'])
        for moltype in moltypes:
            molecule_group = _create_molecule_group(moltype, particle_arrays)
            mol_hierarchy[moltype] = molecule_group

        return mol_hierarchy
