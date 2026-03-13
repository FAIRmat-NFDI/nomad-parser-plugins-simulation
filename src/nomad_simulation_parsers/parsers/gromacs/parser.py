import os
from typing import Any

import numpy as np
from ase.symbols import symbols2numbers
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser.mapping_parser import (
    MappingParser,
    MetainfoParser,
    TextParser,
)
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.force_field import ParticleParametersContainer
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
)
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    MolecularDynamics,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser
from nomad_simulation_parsers.schema_packages import gromacs

from .edr_parser import GromacsEDRParser as GromacsEDRFileParser
from .log_parser import GromacsLogParser as GromacsLogTextParser
from .mdanalysis_parser import GromacsMDAnalysisParser as GromacsMDAnalysisFileParser
from .mdp_parser import GromacsMdpParser as GromacsMDPTextParser
from .xvg_parser import GromacsXvgParser as GromacsXVGTextParser

LOGGER = get_logger(__name__)
ENERGY_UNIT = ureg.kilojoule / ureg.avogadro_number


class GromacsMetainfoParser(MetainfoParser):
    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GromacsThermodynamicsParser(MappingParser):
    _trajectory_steps: list[int]
    _thermodynamic_steps: list[int]

    def __init__(self, **kwargs):
        self._trajectory_steps: list[int] = []
        self._thermodynamic_steps: list[int] = []
        super().__init__(**kwargs)

    _base_calc_unit_map = {
        'Temperature': ureg.kelvin,
        'Volume': ureg.nm**3,
        'Density': ureg.kilogram / ureg.m**3,
        'Pressure (bar)': ureg.bar,
        'Pressure': ureg.bar,
        'Enthalpy': ENERGY_UNIT,
    }
    _energy_label_map = {
        'Potential': 'potential',
        'Kinetic En.': 'kinetic',
        'Total Energy': 'total',
        'pV': 'pressure_volume_work',
    }

    def get_outputs(self, source: dict[str, Any] = {}) -> list[dict[str, Any]]:
        outputs = []
        if not source:
            source = self.data
        times = source.get('Time', [])
        outputs = []
        for n, step in enumerate(self._thermodynamic_steps):
            data = dict(
                step=step,
                time=times[n] * ureg.picosecond if times[n] is not None else None,
            )
            for key in source.keys():
                val = source.get(key)
                if val is None or val[n] is None:
                    continue
                if key in self._energy_label_map:
                    val_n = val[n] * ENERGY_UNIT
                    if key == 'Total Energy':
                        data.setdefault('energy', {})['value'] = val_n
                    else:
                        data.setdefault('energy', {}).setdefault(
                            'contributions', []
                        ).append(dict(value=val_n, label=self._energy_label_map[key]))
                elif key in self._base_calc_unit_map:
                    data[key] = val[n] * self._base_calc_unit_map[key]
            if step in self._trajectory_steps:
                data['system_ref'] = (
                    f'/data/model_system/{self._trajectory_steps.index(step)}'
                )
            outputs.append(data)
        return outputs

    def get_energies(self) -> list[float] | None:
        energies = []
        outputs = self.get_outputs()
        for output in outputs:
            energy = output.get('energy', {}).get('value')
            if energy is None:
                for contribution in output.get('energy', {}).get('contributions', []):
                    if contribution.get('label') == 'potential':
                        energy = contribution.get('value')
                        break
            if energy is not None:
                energies.append(energy)
        if len(energies) != len(self._thermodynamic_steps):
            return None
        return energies


class GromacsLogParser(TextParser, GromacsThermodynamicsParser):
    _trajectory_steps_sampled: list[int]

    def __init__(self, **kwargs):
        self._trajectory_steps_sampled: list[int] = []
        super().__init__(**kwargs)

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self):
        data_object = super().load_file()

        def normalize(input_dict: dict[str, Any]) -> None:
            for key, val in input_dict.items():
                if isinstance(val, dict):
                    normalize(val)
                elif isinstance(val, str):
                    input_dict[key.replace('_', '-')] = val.lower()
                elif isinstance(val, float):
                    if abs(val) == np.inf:
                        input_dict[key] = 'inf' if val > 0 else '-inf'

        input_parameters = data_object.get('input_parameters') or {}
        normalize(input_parameters)
        data_object._results['input_parameters'] = input_parameters

        return data_object

    def get_version(self, source: str) -> str:
        return (source or 'unknown').lstrip('VERSION')

    def get_configurations(self):
        pbc = [
            k in self.data.get('input_parameters', {}).get('pbc', 'xyz') for k in 'xyz'
        ]
        return [dict(pbc=pbc) for _ in self._trajectory_steps_sampled]

    def get_outputs(self) -> list[dict[str, Any]]:
        data = {}
        steps = self.data.get('step', [])
        n_steps = len(self._thermodynamic_steps)
        for n, step in enumerate(steps):
            energies = step.get('energies')
            info = step.get('step_info')
            if energies is None or info is None:
                continue
            step_n = int(info.get('Step'))
            if step_n not in self._thermodynamic_steps:
                continue
            index = self._thermodynamic_steps.index(step_n)
            for key in energies.keys():
                data.setdefault(key, [None] * n_steps)
                data[key][index] = energies.get(key)
            data.setdefault('Time', [None] * n_steps)
            data['Time'][index] = info.get('Time', 1)
        outputs = super().get_outputs(data)
        return outputs

    def get_integrator_type(self, integrator: str) -> str:
        integrator = (integrator or 'md').lower()
        integrator_map = {
            'steep': 'steepest_descent',
            'cg': 'conjugant_gradient',
            'l-bfgs': 'low_memory_broyden_fletcher_goldfarb_shanno',
            'md': 'leap_frog',
            'md-vv': 'velocity_verlet',
            'sd': 'langevin_goga',
            'bd': 'brownian',
        }
        value = integrator_map.get(
            integrator,
            [val for key, val in integrator_map.items() if key in integrator],
        )
        return (
            value
            if not isinstance(value, list)
            else value[0]
            if len(value) != 0
            else None
        )

    def get_coulomb_type(self, coulombtype: str) -> str:
        """Map GROMACS coulombtype to NOMAD schema enum."""
        result = None
        if not coulombtype:
            return result

        coulombtype_lower = coulombtype.lower().replace('_', '-')
        coulomb_map = {
            'cut-off': 'cutoff',
            'cutoff': 'cutoff',
            'ewald': 'ewald',
            'pme': 'particle_mesh_ewald',
            'p3m-ad': 'particle_particle_particle_mesh',
            'reaction-field': 'reaction_field',
            'reaction-field-zero': 'reaction_field',
        }
        result = coulomb_map.get(coulombtype_lower)
        return result

    def get_coordinate_save_frequency(self, input_params: dict[str, Any]) -> int | None:
        """Get coordinate save frequency, preferring compressed output."""
        result = None
        if not input_params:
            return result

        # Prefer nstxout-compressed over nstxout
        nstxout_compressed = input_params.get('nstxout-compressed')
        nstxout = input_params.get('nstxout')

        if nstxout_compressed is not None and nstxout_compressed > 0:
            result = nstxout_compressed
        elif nstxout is not None and nstxout > 0:
            result = nstxout

        return result

    def get_thermodynamic_ensemble(self, input_params: dict[str, Any]) -> str | None:
        """Determine thermodynamic ensemble from thermostat and barostat settings."""
        result = None
        if input_params is None:
            return result

        tcoupl = (input_params.get('tcoupl', '') or 'no').lower()
        pcoupl = (input_params.get('pcoupl', '') or 'no').lower()

        has_thermostat = tcoupl != 'no'
        has_barostat = pcoupl != 'no'

        if has_thermostat and has_barostat:
            result = 'NPT'
        elif has_thermostat:
            result = 'NVT'
        elif has_barostat:
            result = 'NPH'
        else:
            result = 'NVE'

        return result

    def get_thermostat_type(self, tcoupl: str) -> str | None:
        """Map GROMACS tcoupl to NOMAD thermostat_type enum."""
        result = None
        if not tcoupl:
            return result

        tcoupl_lower = tcoupl.lower()
        if tcoupl_lower == 'no':
            return result

        thermostat_map = {
            'berendsen': 'berendsen',
            'nose-hoover': 'nose_hoover',
            'v-rescale': 'velocity_rescaling',
            'andersen': 'andersen',
            'andersen-massive': 'andersen_massive',
        }
        result = thermostat_map.get(tcoupl_lower)
        return result

    def get_reference_temperature(self, input_params: dict[str, Any]) -> float | None:
        """Extract reference temperature from ref_t (handle scalar or array)."""
        result = None
        if not input_params:
            return result

        ref_t = input_params.get('ref-t') or input_params.get('ref_t')
        if ref_t is None:
            return result

        # Handle array (take first value) or scalar
        if isinstance(ref_t, list):
            result = ref_t[0] if len(ref_t) > 0 else None
        else:
            result = ref_t

        return result

    def get_thermostat_coupling_constant(
        self, input_params: dict[str, Any]
    ) -> float | None:
        """Extract thermostat coupling constant from tau_t (handle scalar or array)."""
        result = None
        if not input_params:
            return result

        tau_t = input_params.get('tau-t') or input_params.get('tau_t')
        if tau_t is None:
            return result

        # Handle array (take first value) or scalar
        if isinstance(tau_t, list):
            result = tau_t[0] if len(tau_t) > 0 else None
        else:
            result = tau_t

        return result

    def get_barostat_type(self, pcoupl: str) -> str | None:
        """Map GROMACS pcoupl to NOMAD barostat_type enum."""
        result = None
        if not pcoupl:
            return result

        pcoupl_lower = pcoupl.lower()
        if pcoupl_lower == 'no':
            return result

        barostat_map = {
            'berendsen': 'berendsen',
            'parrinello-rahman': 'parrinello_rahman',
            'mttk': 'mttk',
            'c-rescale': 'c_rescale',
        }
        result = barostat_map.get(pcoupl_lower)
        return result

    def get_barostat_coupling_type(self, pcoupltype: str) -> str | None:
        """Map GROMACS pcoupltype to NOMAD coupling_type enum."""
        result = None
        if not pcoupltype:
            return result

        pcoupltype_lower = pcoupltype.lower()
        coupling_map = {
            'isotropic': 'isotropic',
            'semiisotropic': 'semi_isotropic',
            'anisotropic': 'anisotropic',
            'surface-tension': 'surface_tension',
        }
        result = coupling_map.get(pcoupltype_lower)
        return result

    def get_reference_pressure(self, input_params: dict[str, Any]) -> float | None:
        """Extract reference pressure from ref_p (handle scalar or matrix)."""
        result = None
        if not input_params:
            return result

        ref_p = input_params.get('ref-p') or input_params.get('ref_p')
        if ref_p is None:
            return result

        # Handle matrix/array (take first value) or scalar
        if isinstance(ref_p, list):
            if isinstance(ref_p[0], list):
                # Matrix: take [0][0]
                result = ref_p[0][0] if len(ref_p) > 0 and len(ref_p[0]) > 0 else None
            else:
                # Array: take first value
                result = ref_p[0] if len(ref_p) > 0 else None
        else:
            result = ref_p

        return result

    def get_barostat_coupling_constant(
        self, input_params: dict[str, Any]
    ) -> float | None:
        """Extract barostat coupling constant from tau_p."""
        result = None
        if not input_params:
            return result

        tau_p = input_params.get('tau-p') or input_params.get('tau_p')
        if tau_p is not None:
            result = tau_p

        return result

    def get_compressibility(self, input_params: dict[str, Any]) -> float | None:
        """Extract compressibility (handle scalar or matrix)."""
        result = None
        if not input_params:
            return result

        compressibility = input_params.get('compressibility')
        if compressibility is None:
            return result

        # Handle matrix/array (take first value) or scalar
        if isinstance(compressibility, list):
            if len(compressibility) == 0:
                result = None
            elif isinstance(compressibility[0], list):
                # Matrix: take [0][0]
                result = compressibility[0][0] if len(compressibility[0]) > 0 else None
            else:
                # Array: take first value
                result = compressibility[0]
        else:
            result = compressibility

        return result


class GromacsMDPParser(TextParser):
    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self):
        data_object = super().load_file()

        def check_input_parameters_dict_recursive(input_dict: dict, key: str):
            if key in input_dict:
                return True
            for _, v in input_dict.items():
                if isinstance(v, dict):
                    if check_input_parameters_dict_recursive(v, key):
                        return True
            return False

        input_parameters = data_object.get('input_parameters') or {}
        # parameters that are unique to the mdp file
        input_parameters['mdp_unique_params'] = {}
        for key, param in input_parameters.items():
            new_key = key.replace('_', '-')
            if not check_input_parameters_dict_recursive(input_parameters, new_key):
                input_parameters['mdp_unique_params'][new_key] = (
                    param.lower() if isinstance(param, str) else param
                )

        if data_object._results is None:
            data_object._results = {}
        data_object._results['input_parameters'] = input_parameters
        return data_object


class GromacsEDRParser(GromacsThermodynamicsParser):
    # Fallback parser
    edr_parser: GromacsEDRFileParser = None

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs) -> dict[str | int, Any]:
        if self.data_object is not None:
            for key in self.data_object.keys():
                self.data_object.parse(key)
        return self.data_object._results

    def from_dict(self, dct: dict[str, Any]):
        raise NotImplementedError

    def load_file(self) -> Any:
        if self.filepath:
            self.edr_parser.mainfile = self.filepath
        return self.edr_parser


class GromacsXVGParser(TextParser):
    input_parameters = {}

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_results(self) -> dict[str, Any]:
        title = self.data.get('title', '')
        # TODO incorporate x and y axis labels into the checks
        if not (r'dH/d\xl\f{}' in title and r'\xD\f{}H' in title):
            return {}
        results = {}

        free_energy = results.setdefault('free_energy_calculations', {})
        columns = self.data.get('column_vals')
        free_energy['n_frames'] = len(columns)
        free_energy['value_unit'] = str(ENERGY_UNIT.units)
        # TODO get n_states from input_parameters
        free_energy['n_states'] = self.input_parameters.get('n_states')
        xaxis = self.data.get('xaxis', '').lower()
        # The expected columns of the xvg file are:
        # Total Energy
        # dH/dlambda current lambda
        # Delta H between each lambda and current lambda (n_lambda columns)
        # PV Energy
        if 'time' in xaxis and columns[:, 3:-1].shape[1] == free_energy['n_states']:
            free_energy['times'] = columns[:, 0] * ureg.ps
            columns = columns[:, 1:] * ENERGY_UNIT.magnitude

        return results


class GromacsMDAnalysisParser(MappingParser):
    aux_files: list[str]
    mdanalysis_parser: GromacsMDAnalysisFileParser
    _trajectory_steps_sampled: list[int]
    _trajectory_steps: list[int]
    _thermodynamic_steps: list[int]
    _subsystems_hierarchy: list[dict[str, Any]]

    def __init__(self, **kwargs):
        self.aux_files: list[str] = []
        self.mdanalysis_parser: GromacsMDAnalysisFileParser = None
        self._trajectory_steps_sampled: list[int] = []
        self._trajectory_steps: list[int] = []
        self._thermodynamic_steps: list[int] = []
        self._subsystems_hierarchy: list[dict[str, Any]] = []
        self._hierarchy_returned: bool = False
        super().__init__(**kwargs)

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs) -> dict[str | int, Any]:
        if self.data_object is not None:
            self.data_object.parse()
            result = self.data_object._results

            # Build complete nested subsystem structure and store as parser attribute
            # get_subsystems_from_dict will access this to extract system hierarchy
            self._subsystems_hierarchy = self.get_sub_systems()
            self._hierarchy_returned = False

            return result
        return {}

    def get_subsystems_from_dict(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        """
        Extract subsystems from source dict.

        Called by MappingParser for both top-level and nested ModelSystems.
        For the first top-level call (representative frame), returns the full
        hierarchy from self._subsystems_hierarchy. Subsequent top-level calls
        (other trajectory frames) return empty — topology is frame-independent.
        For nested calls, returns 'sub_systems' from the source dict directly.

        Args:
            source: Source dictionary containing subsystem data
            **kwargs: Additional arguments (unused)

        Returns:
            List of subsystem dicts from this level
        """
        # Root call: first frame gets the full hierarchy; others get nothing.
        # The check 'n_atoms' in source distinguishes a top-level configuration
        # dict (produced by get_configurations) from a nested subsystem dict.
        if 'n_atoms' in source and 'sub_systems' not in source:
            if self._subsystems_hierarchy and not self._hierarchy_returned:
                self._hierarchy_returned = True
                return self._subsystems_hierarchy
            return []

        # Nested call: propagate already-built sub_systems from the hierarchy dict.
        return source.get('sub_systems', []) if isinstance(source, dict) else []

    def from_dict(self, dct: dict[str, Any]):
        raise NotImplementedError

    def load_file(self) -> GromacsMDAnalysisFileParser:
        if self.filepath:
            self.mdanalysis_parser.mainfile = self.filepath
            if self.aux_files:
                self.mdanalysis_parser.auxilliary_files = self.aux_files
        return self.mdanalysis_parser

    def get_particle_labels(self, index: int = 0) -> list[str]:
        labels = self.data_object.get_atom_labels(index)
        try:
            symbols2numbers(labels)
        # TODO: remove blanked `Exception`, hides too many problems.
        except Exception:
            labels = ['CGX'] * len(labels)
        return labels

    def get_particle_parameters(self) -> list[dict[str, Any]]:
        """Return per-particle topology parameters as a list of dicts.

        Each dict contains:
        - 'label' (str): raw force-field site name (e.g. 'OW', 'HW1', 'CA')
        - 'element' (str): chemical element symbol (e.g. 'O', 'H', 'C'),
          absent for coarse-grained or unresolvable particles
        - 'mass' (float, amu): particle mass when available from MDAnalysis
        - 'charge' (float, elementary charge): partial charge when available
        """
        particles_info = {}
        if self.data_object is not None and hasattr(self.data_object, 'get'):
            particles_info = self.data_object.get('atoms_info', {}) or {}
        raw_names = particles_info.get('names', []) or []
        elements = particles_info.get('elements', []) or []
        # Fall back to get_particle_labels when particles_info is not populated
        # (e.g. in tests)
        if not raw_names and not elements:
            fallback = self.get_particle_labels(0)
            raw_names = fallback
            elements = fallback
        n_particles = max(len(raw_names), len(elements))
        result: list[dict[str, Any]] = []
        for i in range(n_particles):
            label = raw_names[i] if i < len(raw_names) else None
            element = elements[i] if i < len(elements) else None
            entry: dict[str, Any] = {'label': label or element or 'CGX'}
            if element and element != 'CGX':
                entry['element'] = element
            result.append(entry)
        universe = getattr(self.data_object, 'universe', None)
        if universe is not None:
            try:
                for i, mass in enumerate(universe.atoms.masses):
                    result[i]['mass'] = float(mass)
            except Exception:
                pass
            try:
                for i, charge in enumerate(universe.atoms.charges):
                    result[i]['charge'] = float(charge)
            except Exception:
                pass
        return result

    def get_particle_parameters_by_type(self) -> list[dict[str, Any]]:
        """Return per-type particle parameters grouped from per-particle data.

        Calls `get_particle_parameters()` and groups by particle-type label. Each unique
        label produces one dict with keys ``particle_type``, and optionally
        ``partial_charge`` (elementary charge) and ``effective_mass`` (amu),
        taken from the first particle of that type encountered.
        """
        per_particle = self.get_particle_parameters()
        seen: dict[str, dict[str, Any]] = {}
        for entry in per_particle:
            atype = entry['label']
            if atype not in seen:
                type_entry: dict[str, Any] = {'particle_type': atype}
                if 'charge' in entry:
                    type_entry['partial_charge'] = entry['charge']
                if 'mass' in entry:
                    type_entry['effective_mass'] = entry['mass']
                seen[atype] = type_entry
        return list(seen.values())

    def get_outputs(self) -> list[dict[str, Any]]:
        outputs = []
        for step in self._thermodynamic_steps:
            if step not in self._trajectory_steps:
                continue
            n = self._trajectory_steps.index(step)
            data = dict(forces=self.data_object.get_forces(n))
            outputs.append(data)
        return outputs

    def get_bond_list(self) -> np.ndarray | None:
        """
        Extract bond list from MDAnalysis interactions.
        Returns numpy array of shape (n_bonds, 2) with atom indices.
        """
        result = None
        interactions = self.data_object.get_interactions()
        if not interactions:
            return result

        # Filter for bond interactions only
        bonds = []
        n_bond_particles: int = 2
        for interaction in interactions:
            if interaction.get('type') == 'bond':
                particle_indices = interaction.get('atom_indices')
                if (
                    particle_indices is not None
                    and len(particle_indices) == n_bond_particles
                ):
                    bonds.append(list(particle_indices))

        if bonds:
            result = np.array(bonds, dtype=int)

        return result

    def get_configurations(self) -> list[dict[str, Any]]:
        # Labels, n_particles and bond_list are frame-independent; compute once.
        labels = self.get_particle_parameters()
        n_particles = self.data_object.get_n_atoms(0)
        bond_list = self.get_bond_list()

        configurations = []
        for n, _ in enumerate(self._trajectory_steps_sampled):
            # get_frame_data performs a single trajectory seek for all per-frame data.
            frame_data = self.data_object.get_frame_data(n)
            config = dict(
                n_atoms=n_particles,
                labels=labels,
                positions=frame_data['positions'],
                velocities=frame_data['velocities'],
                lattice_vectors=frame_data['lattice_vectors'],
            )
            if n == 0 and bond_list is not None:
                config['bond_list'] = bond_list
            configurations.append(config)
        return configurations

    def get_force_field_contributions(self) -> list[dict[str, Any]]:
        """
        Transform MDAnalysis interactions into force field contributions.
        Groups interactions by type and creates one Potential per interaction type.

        CURRENT IMPLEMENTATION:
        ======================
        - Extracts topology (bond/angle/dihedral connectivity) from MDAnalysis
        - Groups by interaction type (e.g., all bonds together)
        - Provides particle_indices and particle_labels
        - Does NOT include force field parameters

        LIMITATION:
        ==========
        Force field parameters are stored separately in TPR files (see
        get_force_field_parameters() in mdanalysis_parser.py). Connecting them
        requires parsing the interaction lists (ilist) that map each bond/angle
        to its parameter set index.

        Example TPR structure (gromacs/src/gromacs/fileio/tpxio.cpp):
        - functypes = [F_BONDS, F_ANGLES, F_LJ, ...]  # parameter type IDs
        - parameters = [[k1, r01], [k2, r02], ...]     # parameter values
        - ilist[F_BONDS] = [[0, 1, 0], [1, 2, 1], ...]  # [atom_i, atom_j, param_idx]

        For full implementation:
        1. Parse ilist sections from TPR (requires custom reader or MDAnalysis
           extension)
        2. Match each interaction to parameter set: params[ilist[type][i][2]]
        3. Populate Potential.parameters with actual force constants
        4. Map GROMACS function types to NOMAD potential classes:
           - F_BONDS (harmonic) -> HarmonicBond
           - F_MORSE -> MorseBond
           - F_ANGLES -> HarmonicAngle
           - etc.

        NOMAD Schema Mapping:
        ====================
        ForceField.contributions should contain:
        - functional_form: 'harmonic_bond', 'morse_bond', etc.
        - particle_indices: [[i, j], [k, l], ...] for each bond/angle
        - particle_labels: [['O', 'H'], ...] (optional)
        - parameters: Subsection with actual values (e.g., HarmonicBond.bond_constant)

        Current output is suitable for topology visualization but lacks
        force field parameter details.
        """

        interactions = self.data_object.get_interactions()
        if not interactions:
            return []

        # Group interactions by type (bonds, angles, dihedrals, impropers)
        grouped = {}
        for interaction in interactions:
            int_type = interaction.get('type', 'unknown')
            if int_type not in grouped:
                grouped[int_type] = []
            grouped[int_type].append(interaction)

        # Create one Potential contribution per interaction type
        contributions = []
        for int_type, int_list in grouped.items():
            # Collect all indices and labels for this type
            all_indices = []
            all_labels = []
            for inter in int_list:
                particle_indices = inter.get('atom_indices')
                particle_labels = inter.get('atom_labels')
                if particle_indices is not None and len(particle_indices) > 0:
                    all_indices.append(list(particle_indices))
                if particle_labels is not None and len(particle_labels) > 0:
                    all_labels.append(list(particle_labels))

            if len(all_indices) == 0:
                continue

            contribution = {
                'functional_form': int_type,
                'particle_indices': all_indices,
            }
            if len(all_labels) > 0:
                contribution['particle_labels'] = all_labels

            contributions.append(contribution)

        return contributions

    def get_sub_systems(
        self, source: dict[str, Any] = None, **kwargs
    ) -> list[dict[str, Any]]:
        particles_groups = self.mdanalysis_parser.parse_molecular_hierarchy()
        if particles_groups is None:
            self.logger.warning(
                'parse_molecular_hierarchy returned None; skipping system hierarchy'
            )
            return []
        return list(particles_groups.values())


class GromacsArchiveWriter(MDParser):
    def __init__(self, **kwargs):
        self._simulation_parser = GromacsMetainfoParser()
        self._simulation_parser.max_nested_level = 10
        self._log_parser = GromacsLogParser(text_parser=GromacsLogTextParser())
        self._mdp_parser = GromacsMDPParser(text_parser=GromacsMDPTextParser())
        self._edr_parser = GromacsEDRParser(edr_parser=GromacsEDRFileParser())
        self._mdanalysis_parser = GromacsMDAnalysisParser(
            mdanalysis_parser=GromacsMDAnalysisFileParser()
        )
        self._xvg_parser = GromacsXVGParser(text_parser=GromacsXVGTextParser())
        self._mdp_ext = 'mdp'
        self._mdp_std_filename = 'mdout'
        self.input_parameters = {}
        super().__init__(**kwargs)

    def get_mdp_file(self):
        """
        Tries to find the mdp input parameters (ext = mdp) that match the mainfile
        calculation.
        Priority is as follows:
            1. output mdp file containing both the matching mainfile name and the
            standard gromacs name `mdout`
            2. file containing the standard gromacs name `mdout`
            3. input mdp file matching the mainfile name (as usual)
            4. any `.mdp` file within the directory (as usual)
        """
        files = [d for d in self._gromacs_files if d.endswith(self._mdp_ext)]

        if len(files) == 0:
            return ''

        if len(files) == 1:
            return os.path.join(self._maindir, files[0])

        for f in files:
            filename = f.rsplit('.', 1)[0]
            if self._basename in filename and self._mdp_std_filename in filename:
                return os.path.join(self._maindir, f)

        for f in files:
            filename = f.rsplit('.', 1)[0]
            if self._mdp_std_filename in filename:
                return os.path.join(self._maindir, f)

        return self.get_gromacs_file(self._mdp_ext)

    def get_gromacs_file(self, ext: str) -> str:
        files = [d for d in self._gromacs_files if d.endswith(ext)]

        if len(files) == 0:
            return ''

        if len(files) == 1:
            return os.path.join(self._maindir, files[0])

        basenames = [f.rsplit('.', 1)[0] for f in files]
        # we assume that the file has the same basename as the log file e.g.
        # out.log would correspond to out.tpr and out.trr and out.edr
        for n, basename in enumerate(basenames):
            if basename == self._basename:
                return os.path.join(self._maindir, files[n])

        for n, basename in enumerate(basenames):
            if basename.startswith(self._basename):
                return os.path.join(self._maindir, files[n])

        # if the files are all named differently, we guess that the one that does not
        # share the same basename would be file we are interested in
        # e.g. out.log someout.log out.tpr out.trr another.tpr file.trr
        # we guess that the out.* files belong together and the rest that does not share
        # a basename would be grouped together
        counts = []
        all_basenames = [f.rsplit('.', 1)[0] for f in self._gromacs_files]
        for n, basename in enumerate(basenames):
            count = 0
            for ref_basename in all_basenames:
                if basename == ref_basename:
                    count += 1
            if count == 1:
                return os.path.join(self._maindir, files[n])
            counts.append(count)

        return os.path.join(self._maindir, files[counts.index(min(counts))])

    def _parse_workflow_section(self):
        integrator = self.input_parameters.get('integrator', 'md').lower()
        workflow2 = None
        if integrator in ['l-bfgs', 'cg', 'steep']:
            workflow2 = GeometryOptimization()
        else:
            workflow2 = MolecularDynamics()

        if workflow2 is None:
            return
        self._simulation_parser.data_object = workflow2

        # parse main log file
        self._simulation_parser.annotation_key = gromacs.LOG_KEY
        self._log_parser.convert(self._simulation_parser)

        # parse edr file
        self._simulation_parser.annotation_key = gromacs.EDR_KEY
        self._edr_parser.convert(self._simulation_parser)

        self.archive.workflow2 = workflow2

    def _parse_data_section(self):
        self.archive.data = Simulation(program=Program(name='GROMACS'))
        self._simulation_parser.data_object = self.archive.data

        # parse main log file
        self._simulation_parser.annotation_key = gromacs.LOG_KEY
        self._log_parser.convert(self._simulation_parser)

        # parse edr file
        self._simulation_parser.annotation_key = gromacs.EDR_KEY
        self._edr_parser.convert(self._simulation_parser)

        # parse mdanalysis trajectory files
        self._simulation_parser.annotation_key = gromacs.TPR_KEY
        self._mdanalysis_parser.convert(self._simulation_parser)

        # parse per-type particle parameters (charges, masses, particle types) from
        # TPR into a fresh ParticleParametersContainer, then append it to
        # ForceField.numerical_settings. Using a dedicated PARTICLE_PARAM_KEY and a
        # separate data_object avoids the index-0 collision with ForceCalculations
        # written by the LOG pass.
        if self.archive.data.model_method:
            ppc = ParticleParametersContainer()
            self._simulation_parser.data_object = ppc
            self._simulation_parser.annotation_key = gromacs.PARTICLE_PARAM_KEY
            self._mdanalysis_parser.convert(self._simulation_parser)
            if ppc.particle_parameters:
                self.archive.data.model_method[-1].numerical_settings.append(ppc)

    def write_to_archive(self):
        # intitialize variables
        self._maindir = os.path.dirname(self.mainfile)
        self._gromacs_files = os.listdir(self._maindir)
        self._basename = os.path.basename(self.mainfile).rsplit('.', 1)[0]

        # set up source parsers
        self._log_parser.filepath = self.mainfile
        self._edr_parser.filepath = self.get_gromacs_file('edr')
        self._mdanalysis_parser.filepath = self.get_gromacs_file('tpr')
        # TODO include input parameters read from mdp parser
        self._mdp_parser.filepath = self.get_mdp_file()
        self._xvg_parser.filepath = self.get_gromacs_file('xvg')
        # determine auxiliary file order: trr, xtc
        for ext in ['trr', 'xtc']:
            aux_file = self.get_gromacs_file(ext)
            if not aux_file:
                continue
            self._mdanalysis_parser.aux_files = [aux_file]
            # check if positions are parsed
            positions = self._mdanalysis_parser.data_object.get_positions(0)
            if positions is not None:
                break

        # build input parameters from log and mdp
        self.input_parameters = {
            **self._log_parser.data_object.get('input_parameters', {}),
            **self._mdp_parser.data_object.get('input_parameters', {}),
        }

        # determine sampled trajectory steps
        n_frames = self._mdanalysis_parser.data_object.get('n_frames', 0)
        traj_sampling_rate = self.input_parameters.get('nstxout', 1)
        n_particles_per_frame = self._mdanalysis_parser.data_object.get_n_atoms(0)
        self.n_particles = [n_particles_per_frame] * n_frames
        traj_steps = [n * traj_sampling_rate for n in range(n_frames)]
        self.trajectory_steps = traj_steps

        # determine sampled thermodynamic steps
        calculation_times = self._edr_parser.data.get('Time', [])
        time_step = self.input_parameters.get('dt')
        if time_step is None and len(calculation_times) > 1:
            time_step = calculation_times[1] - calculation_times[0]
        self.thermodynamics_steps = [
            int(time / time_step if time_step else 1) for time in calculation_times
        ]

        # pass sampled steps to parsers
        self._log_parser._trajectory_steps_sampled = self.trajectory_steps
        self._log_parser._trajectory_steps = traj_steps
        self._log_parser._thermodynamic_steps = self.thermodynamics_steps
        self._mdanalysis_parser._trajectory_steps_sampled = self.trajectory_steps
        self._mdanalysis_parser._trajectory_steps = traj_steps
        self._mdanalysis_parser._thermodynamic_steps = self.thermodynamics_steps
        self._edr_parser._thermodynamic_steps = self.thermodynamics_steps
        self._edr_parser._trajectory_steps = traj_steps

        self._parse_data_section()

        self._parse_workflow_section()

        for parser in [
            self._simulation_parser,
            self._log_parser,
            self._mdanalysis_parser,
            self._edr_parser,
            self._mdp_parser,
            self._xvg_parser,
        ]:
            parser.close()


class GromacsParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = GromacsArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = None,
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
