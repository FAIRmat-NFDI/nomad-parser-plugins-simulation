import operator
import re
from collections.abc import Callable
from typing import Any

import MDAnalysis
from MDAnalysis.topology.tpr import setting as tpr_setting
from MDAnalysis.topology.tpr import utils as tpr_utils

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
