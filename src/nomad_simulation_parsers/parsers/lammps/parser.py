import os
import re
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad_simulations.schema_packages.force_field import (
    AnglePotential,
    BondPotential,
    DihedralPotential,
    ForceCalculations,
    ForceField,
    ImproperDihedralPotential,
    ParticleParameters,
    ParticleParametersContainer,
)
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem
from nomad_simulations.schema_packages.workflow.general import (
    ForceConvergenceTarget,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
    GeometryOptimizationModel,
    GeometryOptimizationResults,
)
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
from nomad_simulation_parsers.parsers.utils.constants import (
    CHEMICAL_SYMBOLS,
    REFERENCE_MASSES,
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
        self._md_method: MolecularDynamicsMethod | None = None

    _magic_two: int = 2
    _magic_three: int = 3
    _magic_four: int = 4
    _magic_five: int = 5

    _energy_contribution_keys: dict[str, str] = {
        'KinEng': 'kinetic',
        'PotEng': 'potential',
        'E_pair': 'pair',
        'E_bond': 'bond',
        'E_angle': 'angle',
        'E_dihed': 'dihedral',
        'E_impro': 'improper',
        'E_vdwl': 'van_der_waals',
        'E_coul': 'coulomb',
        'E_kspce': 'kspace',
        'E_long': 'long_range',
        'E_tail': 'tail_correction',
    }

    _min_style_map: dict[str, str] = {
        'cg': 'polak_ribiere_conjugant_gradient',
        'sd': 'steepest_descent',
        'hftn': 'hessian_free_truncated_newton',
        'quickmin': 'damped_dynamics',
        'fire': 'damped_dynamics',
        'spin': 'damped_dynamics',
    }

    _kspace_coulomb_prefix_map: dict[str, str] = {
        'ewald': 'ewald',
        'pppm': 'particle_mesh_ewald',
        'msm': 'multilevel_summation',
    }

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

    def _masses_from_data_file(self) -> dict[int, float]:
        """Parse Masses section from data file → {type_id: mass}."""
        type_mass: dict[int, float] = {}
        masses_data = self._data_parser.get('Masses', None)
        if not masses_data or not isinstance(masses_data, list):
            return type_mass
        raw = masses_data[0][1] if isinstance(masses_data[0], tuple) else masses_data
        if raw is None or len(raw) == 0:
            return type_mass
        for row in raw:
            try:
                type_mass[int(row[0])] = float(row[1])
            except (ValueError, TypeError, IndexError):
                continue
        return type_mass

    def _apply_mass_commands(self, type_mass: dict[int, float]) -> None:
        """Apply `mass T M` input-script commands to the type→mass map."""
        mass_cmds = self._log_parser.get('mass')
        if mass_cmds is None:
            return
        if isinstance(mass_cmds, list) and mass_cmds and isinstance(mass_cmds[0], str):
            mass_cmds = [mass_cmds]
        for cmd in mass_cmds:
            try:
                if len(cmd) >= self._magic_two:
                    type_mass[int(float(cmd[0]))] = float(cmd[1])
            except (ValueError, TypeError, IndexError):
                continue

    @staticmethod
    def _parse_set_mass_cmd(
        cmd: list[str],
    ) -> tuple[str | None, str | None, float | None]:
        """Return (style, target, mass) from a parsed `set` token list.

        Scans keyword-value pairs after the style and target tokens.
        Returns (None, None, None) if the command does not set mass.
        """
        _MIN_LEN = 3
        if not isinstance(cmd, list) or len(cmd) < _MIN_LEN:
            return None, None, None
        style = cmd[0]
        target = cmd[1]
        i = _MIN_LEN - 1
        while i < len(cmd) - 1:
            if cmd[i] == 'mass':
                try:
                    return style, target, float(cmd[i + 1])
                except (ValueError, IndexError):
                    return style, target, None
            i += 2
        return style, target, None

    def _apply_set_type_masses(self, type_mass: dict[int, float]) -> None:
        """Apply `set type T mass M` commands; warn on per-atom overrides."""
        set_cmds = self._log_parser.get('set')
        if set_cmds is None:
            return
        if isinstance(set_cmds, list) and set_cmds and isinstance(set_cmds[0], str):
            set_cmds = [set_cmds]
        for cmd in set_cmds:
            style, target, mass_val = self._parse_set_mass_cmd(cmd)
            if style is None or mass_val is None:
                continue
            if style == 'type':
                try:
                    type_mass[int(float(target))] = mass_val
                except (ValueError, TypeError):
                    pass
            elif style in ('atom', 'group', 'region'):
                # TODO: implement virtual-type mass overrides for `set atom`
                # For `set atom ID mass M`:
                #   1. Build atom_id → type_id from the data file Atoms
                #      section (_atom_id_to_type helper).
                #   2. Allocate a synthetic type_id (e.g. max_real_type + n)
                #      that is unique within this simulation.
                #   3. Add {synthetic_type_id: mass_val} to type_mass.
                #   4. Store {atom_id: synthetic_type_id} in a new instance
                #      dict (e.g. self._atom_type_overrides) so it can be
                #      passed to TrajParser via _set_parser_masses.
                # `set group` / `set region` require parsing group membership,
                # which is only feasible for simple id-range groups; leave
                # those with the existing warning for now.
                self.logger.warning(
                    'Per-atom mass override not applied to per-type map',
                    set_style=style,
                    target=target,
                )

    def _masses_from_eam_bop(self) -> None:
        """Stub for EAM/BOP potential-file mass extraction.

        EAM and BOP potential files may embed masses that take precedence
        over data/input-script values. Parsing those files is not yet
        implemented; a warning is emitted and the caller falls back to
        whatever masses were already collected.
        """
        pair_style_raw = self._log_parser.get('pair_style')
        if pair_style_raw is None:
            return
        ps = pair_style_raw[0] if isinstance(pair_style_raw, list) else pair_style_raw
        if isinstance(ps, str) and ps.lower() in ('eam', 'bop'):
            self.logger.warning(
                'EAM/BOP potential files may embed masses; '
                'potential-file mass override not implemented, '
                'using data/input-script masses as fallback'
            )

    def _extract_masses(self) -> np.ndarray | None:
        """Build a per-type mass array with sequential override semantics.

        Priority order (later wins):
        1. Data file ``Masses`` section
        2. ``mass T M`` commands in the input script
        3. ``set type T mass M`` commands in the input script

        ``set atom/group/region mass`` overrides are noted via a warning
        but are not applied to the per-type map. EAM/BOP potential-file
        masses are not implemented; a warning is emitted if detected.

        Returns an ``(n_types, 2)`` array ``[[type_id, mass], ...]``,
        or ``None`` if no mass data was found.
        """
        type_mass = self._masses_from_data_file()
        self._apply_mass_commands(type_mass)
        self._apply_set_type_masses(type_mass)
        self._masses_from_eam_bop()
        if not type_mass:
            return None
        return np.array([[tid, m] for tid, m in sorted(type_mass.items())])

    def _extract_per_type_charges(self) -> dict[int, float]:
        """Return a {type_id: charge} map from the Atoms section.

        Only populated for atom_styles 'full' and 'charge'; empty otherwise.
        """
        atom_style_raw = self._log_parser.get('atom_style')
        atom_style = (
            atom_style_raw[0]
            if isinstance(atom_style_raw, list) and atom_style_raw
            else (atom_style_raw or '')
        )
        per_type_charge: dict[int, float] = {}
        if atom_style not in ('full', 'charge'):
            return per_type_charge
        type_col, charge_col = (2, 3) if atom_style == 'full' else (1, 2)
        atoms_data = self._data_parser.get('Atoms', None)
        if atoms_data is None or len(atoms_data) == 0:
            return per_type_charge
        atoms_arr = atoms_data[0][1] if isinstance(atoms_data[0], tuple) else None
        if (
            atoms_arr is not None
            and atoms_arr.ndim == self._magic_two
            and atoms_arr.shape[1] > charge_col
        ):
            for row in atoms_arr:
                t = int(row[type_col])
                if t not in per_type_charge:
                    per_type_charge[t] = float(row[charge_col])
        return per_type_charge

    def _get_particle_parameters_by_type(
        self,
    ) -> list[dict[str, Any]]:
        """Return per-type particle parameters.

        One dict per unique particle label with keys ``particle_type``
        (matching ``ps.label`` in either ``AtomsState`` or
        ``CGBeadState``), ``effective_mass`` (pint Quantity), and
        optionally ``partial_charge`` (pint Quantity).

        The label mapping is taken from ``TrajParser.chemical_symbols``
        when available so that ``particle_type`` is always consistent
        with what ``get_atom_labels`` returns, regardless of whether
        the particles are atomic or coarse-grained.
        """
        masses_array = self._extract_masses()
        if masses_array is None or len(masses_array) == 0:
            return []

        mass_unit = self._log_parser.units.get('mass', ureg.amu)
        charge_unit = self._log_parser.units.get('charge', ureg.elementary_charge)

        # Use the label mapping already computed by TrajParser so that
        # particle_type agrees exactly with ps.label.  TrajParser keys
        # are np.float64 values; probe with both the raw array value and
        # a plain float to handle both representations.
        # chemical_symbols is only defined on TrajParser, not MDAnalysisParser;
        # access it directly to avoid a spurious "not found in parsers" warning.
        chem_sym = None
        for _p in self.traj_parsers._parsers:
            if isinstance(_p, TrajParser):
                chem_sym = _p.chemical_symbols
                if chem_sym is not None:
                    break

        type_to_label: dict[int, str] = {}
        for row in masses_array:
            type_id = int(row[0])
            mass_val = float(row[1])
            label = None
            if chem_sym is not None:
                label = chem_sym.get(row[0]) or chem_sym.get(float(type_id))
            if label is None:
                symbol_idx = int(np.argmin(np.abs(REFERENCE_MASSES - mass_val)))
                label = CHEMICAL_SYMBOLS[symbol_idx]
            if label == 'X':
                label = 'CGX'
            type_to_label[type_id] = label

        per_type_charge = self._extract_per_type_charges()

        # One ParticleParameters entry per unique label.  Multiple type
        # IDs that resolve to the same label share one entry because
        # ps.label is also the same for all of them.
        seen: dict[str, dict[str, Any]] = {}
        for row in masses_array:
            type_id = int(row[0])
            mass_val = float(row[1])
            label = type_to_label[type_id]
            if label not in seen:
                entry: dict[str, Any] = {
                    'particle_type': label,
                    'effective_mass': mass_val * mass_unit,
                }
                if type_id in per_type_charge:
                    entry['partial_charge'] = per_type_charge[type_id] * charge_unit
                seen[label] = entry
        return list(seen.values())

    def _set_parser_masses(self, masses: np.ndarray | None) -> None:
        """Set masses on trajectory parsers."""
        for parser in self.traj_parsers._parsers:
            if isinstance(parser, TrajParser):
                parser.masses = masses
                # TODO: also pass self._atom_type_overrides here once
                # virtual-type support is implemented in
                # _apply_set_type_masses.  TrajParser will need a new
                # `atom_type_overrides` attribute (dict[int, int]) that
                # get_atom_labels checks before the chemical_symbols
                # lookup, overriding the type_id for specific atom IDs.

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

    def _build_md_method(self) -> MolecularDynamicsMethod:
        """Build and return a MolecularDynamicsMethod from log file data."""
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

        return method

    def _extract_pair_cutoffs(self) -> tuple[float | None, float | None]:
        """Return (vdw_cutoff, coulomb_cutoff) in raw distance units.

        Takes the last trailing numeric token from pair_style args as the
        outer cutoff. When the style name contains 'coul/long' or 'coul/msm',
        that cutoff applies to both VDW and Coulomb real-space.
        """
        pair_raw = self._log_parser.get('pair_style')
        if pair_raw is None:
            return None, None
        args = pair_raw[0] if isinstance(pair_raw, list) else pair_raw
        if not isinstance(args, list):
            args = [args] if args else []
        style = args[0].lower() if args else ''
        cutoffs = []
        for tok in args[1:]:
            try:
                cutoffs.append(float(tok))
            except (ValueError, TypeError):
                break
        if not cutoffs:
            return None, None
        vdw_cut = cutoffs[-1]
        has_long_coul = 'coul/long' in style or 'coul/msm' in style
        coul_cut = cutoffs[-1] if has_long_coul else None
        if 'coul/cut' in style and len(cutoffs) >= self._magic_two:
            coul_cut = cutoffs[-1]
        return vdw_cut, coul_cut

    def _extract_kspace_coulomb_type(self) -> str | None:
        """Map kspace_style keyword to ForceCalculations.coulomb_type enum."""
        kspace_raw = self._log_parser.get('kspace_style')
        if kspace_raw is None:
            return None
        args = kspace_raw[0] if isinstance(kspace_raw, list) else kspace_raw
        style = (args[0] if isinstance(args, list) else str(args)).lower()
        prefix = style.split('/')[0]
        return self._kspace_coulomb_prefix_map.get(prefix)

    def _extract_neighbor_skin(self) -> float | None:
        """Return the neighbor list skin distance in raw distance units.

        The ``neighbor skin style`` command sets an extra distance added on
        top of the pair force cutoff. LAMMPS stores all atom pairs within
        ``pair_force_cutoff + skin`` in the neighbor list. The skin lets the
        list go stale for several timesteps before a rebuild is needed.
        The caller adds this value to the pair cutoff to obtain the full
        ``neighbor_update_cutoff`` stored in the schema.
        """
        neighbor_raw = self._log_parser.get('neighbor')
        if neighbor_raw is None:
            return None
        args = neighbor_raw[0] if isinstance(neighbor_raw, list) else neighbor_raw
        tok = args[0] if isinstance(args, list) else args
        try:
            return float(tok)
        except (ValueError, TypeError):
            return None

    def _extract_neigh_modify_frequency(self) -> int | None:
        """Return the 'every N' value from neigh_modify."""
        nm_raw = self._log_parser.get('neigh_modify')
        if nm_raw is None:
            return None
        if isinstance(nm_raw, list) and nm_raw and isinstance(nm_raw[0], str):
            nm_raw = [nm_raw]
        cmd = nm_raw[-1] if isinstance(nm_raw, list) else nm_raw
        if not isinstance(cmd, list):
            return None
        result = None
        for i in range(len(cmd) - 1):
            if cmd[i] == 'every':
                try:
                    result = int(cmd[i + 1])
                except (ValueError, TypeError):
                    pass
                break
        return result

    def _build_force_calculations(self) -> ForceCalculations | None:
        """Build ForceCalculations from pair_style, kspace_style, neighbor list."""
        vdw_cut, coul_cut = self._extract_pair_cutoffs()
        kspace_type = self._extract_kspace_coulomb_type()
        skin = self._extract_neighbor_skin()
        neigh_freq = self._extract_neigh_modify_frequency()
        if all(v is None for v in (vdw_cut, coul_cut, kspace_type, skin, neigh_freq)):
            return None
        fc = ForceCalculations()
        if vdw_cut is not None:
            fc.vdw_cutoff = self.apply_unit(vdw_cut, 'distance')
        if coul_cut is not None:
            fc.coulomb_cutoff = self.apply_unit(coul_cut, 'distance')
        if kspace_type is not None:
            fc.coulomb_type = kspace_type
        outer = vdw_cut if vdw_cut is not None else (coul_cut or 0.0)
        if skin is not None:
            fc.neighbor_update_cutoff = self.apply_unit(outer + skin, 'distance')
        if neigh_freq is not None:
            fc.neighbor_update_frequency = neigh_freq
        return fc

    def _get_style_functional_form(self, style_key: str) -> str | None:
        """Return the functional form string from a LAMMPS style command.

        Single-keyword styles return the keyword directly (e.g. 'harmonic').
        Multi-keyword styles (hybrid ...) are joined with '/' so the full
        specification is preserved (e.g. 'hybrid/harmonic/cosine').
        """
        raw = self._log_parser.get(style_key)
        if raw is None:
            return None
        val = raw[-1] if isinstance(raw, list) else raw
        if isinstance(val, list):
            return '/'.join(val)
        return str(val)

    def _build_force_field_contributions(self) -> list:
        """Build Potential instances from data file topology sections.

        Reads Bonds, Angles, Dihedrals, and Impropers sections and maps each
        to a typed Potential subclass with particle_indices (0-based) and
        functional_form from the corresponding style command in the log.
        """
        _sections: list[tuple] = [
            ('Bonds', BondPotential, 'bond_style', 2, 4),
            ('Angles', AnglePotential, 'angle_style', 2, 5),
            ('Dihedrals', DihedralPotential, 'dihedral_style', 2, 6),
            ('Impropers', ImproperDihedralPotential, 'improper_style', 2, 6),
        ]
        contributions = []
        for section, cls, style_key, col_start, col_end in _sections:
            data = self._data_parser.get(section, None)
            if data is None or len(data) == 0:
                continue
            array = data[0][1] if isinstance(data[0], tuple) else None
            if (
                array is None
                or array.ndim != self._magic_two
                or array.shape[1] < col_end
            ):
                continue
            pot = cls()
            functional_form = self._get_style_functional_form(style_key)
            if functional_form is not None:
                pot.functional_form = functional_form
            pot.particle_indices = array[:, col_start:col_end].astype(np.int32) - 1
            contributions.append(pot)
        return contributions

    def _build_force_field(self) -> ForceField | None:
        """Build a ForceField with topology contributions and particle parameters.

        Returns None when no data is available from any source.
        """
        pp_by_type = self._get_particle_parameters_by_type()
        force_calcs = self._build_force_calculations()
        contributions = self._build_force_field_contributions()

        if not pp_by_type and force_calcs is None and not contributions:
            return None

        force_field = ForceField()
        ppc = ParticleParametersContainer()
        for entry in pp_by_type:
            pp = ParticleParameters()
            pp.particle_type = entry['particle_type']
            if 'effective_mass' in entry:
                pp.effective_mass = entry['effective_mass']
            if 'partial_charge' in entry:
                pp.partial_charge = entry['partial_charge']
            ppc.particle_parameters.append(pp)

        force_field.numerical_settings.append(ppc)

        if force_calcs is not None:
            force_field.numerical_settings.append(force_calcs)

        for pot in contributions:
            force_field.contributions.append(pot)

        return force_field

    def parse_method(self, simulation: Simulation) -> None:
        """
        Parse method information from data file and log file.

        Populates simulation.model_method with a ForceField (particle
        parameters) and stores MolecularDynamicsMethod in self._md_method
        for later attachment to the workflow section.

        TODO: Still to migrate from legacy parser
        - Parse force calculation parameters (pair_style, kspace_style,
          neighbor) into ForceField.numerical_settings as ForceCalculations.
        """
        if self.traj_parsers[0].mainfile is None or self._data_parser.mainfile is None:
            return

        if self.traj_parsers.eval('n_frames') is None:
            return

        self._md_method = self._build_md_method()

        masses = self._extract_masses()
        self._set_parser_masses(masses)

        force_field = self._build_force_field()
        if force_field is not None:
            simulation.model_method.append(force_field)

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
                # Convert from 1-based (LAMMPS data file) to 0-based indexing (NOMAD)
                # Extract columns 2:4 (atom IDs), convert to int, subtract 1
                bond_array = bonds[0][1][:, 2:4].astype(int) - 1
                self._bond_list = bond_array

        for step in self.trajectory_steps:
            traj_n = self._trajectory_steps.index(step)

            # Extract and apply units to trajectory data
            lattice_vectors = _get_quantity_with_units(
                traj_n, 'lattice_vectors', 'distance'
            )
            if lattice_vectors is None or np.allclose(
                lattice_vectors.magnitude, 0, atol=1e-10
            ):
                lattice_vectors = self._data_parser.get_lattice_vectors()
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
        self,
        molecule: int,
        i_molecule: int,
        particle_arrays: dict,
        moltype_name: str = '',
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
            name=moltype_name if moltype_name else molecule,
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
                mol = self._create_molecule(
                    molecule, i_molecule, particle_arrays, moltype_name=str(moltype)
                )
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

        Extracts per-step thermo data (energy, temperature) and creates a
        TrajectoryOutputs section for each step. model_system_ref is linked
        automatically by parse_output_step when the step number matches a
        sampled trajectory step.
        """
        thermo_data = self._log_parser.get_thermodynamic_data()
        if thermo_data is None:
            return

        steps = thermo_data.get('Step')
        if steps is None or len(steps) == 0:
            return

        self.thermodynamics_steps = [int(s) for s in steps]

        timestep = self._extract_timestep()

        for n, step in enumerate(self.thermodynamics_steps):
            data: dict[str, Any] = {'step': step}

            if timestep is not None:
                data['time'] = float(step) * timestep

            energy_data: dict[str, Any] = {}
            tot_eng = thermo_data.get('TotEng')
            if tot_eng is not None and n < len(tot_eng):
                energy_data['value'] = tot_eng[n]

            contributions = []
            for lammps_key, label in self._energy_contribution_keys.items():
                vals = thermo_data.get(lammps_key)
                if vals is not None and n < len(vals):
                    contributions.append({'value': vals[n], 'label': label})
            if contributions:
                energy_data['contributions'] = contributions
            if energy_data:
                data['total_energies'] = energy_data

            temp = thermo_data.get('Temp')
            if temp is not None and n < len(temp):
                data['temperatures'] = {'value': temp[n]}

            self.parse_output_step(data, simulation)

    @staticmethod
    def _parse_minimization_stats(
        stats: str,
    ) -> dict[str, float | None]:
        result: dict[str, float | None] = {
            'final_energy_difference': None,
            'final_force_maximum': None,
        }
        energy_match = re.search(
            r'Energy initial, next-to-last, final\s*=\s*'
            r'[\d\.\+\-eE]+\s+([\d\.\+\-eE]+)\s+([\d\.\+\-eE]+)',
            stats,
        )
        if energy_match:
            nml = float(energy_match.group(1))
            final = float(energy_match.group(2))
            result['final_energy_difference'] = abs(nml - final)
        force_match = re.search(
            r'Force max component initial, final\s*='
            r'\s*[\d\.\+\-eE]+\s+([\d\.\+\-eE]+)',
            stats,
        )
        if force_match:
            result['final_force_maximum'] = float(force_match.group(1))
        return result

    def _parse_geometry_optimization_workflow(self) -> None:
        model = GeometryOptimizationModel()
        model.optimization_type = 'atomic'

        min_style = self._log_parser.get('min_style')
        if min_style is not None and len(min_style) > 0:
            style = min_style[0]
            model.optimization_method = self._min_style_map.get(style, style)

        minimize_cmd = self._log_parser.get('minimize') or self._log_parser.get(
            'minimize/kk'
        )
        if minimize_cmd is not None and len(minimize_cmd) > 0:
            args = minimize_cmd[0]
            if isinstance(args, list) and len(args) >= self._magic_three:
                try:
                    ftol = float(args[1])
                    maxiter = int(float(args[2]))
                    model.n_steps_maximum = maxiter
                    if ftol > 0.0:
                        force_unit = self._log_parser.units.get('force', 1)
                        fc_target = ForceConvergenceTarget()
                        fc_target.threshold = ftol * force_unit
                        fc_target.threshold_type = 'maximum'
                        model.convergence_targets.append(fc_target)
                except (ValueError, IndexError):
                    pass

        results = GeometryOptimizationResults()
        min_stats = self._log_parser.get('minimization_stats')
        if min_stats is not None:
            parsed = self._parse_minimization_stats(min_stats)
            energy_unit = self._log_parser.units.get('energy', 1)
            force_unit = self._log_parser.units.get('force', 1)
            de = parsed.get('final_energy_difference')
            if de is not None:
                results.final_energy_difference = de * energy_unit
            fmax = parsed.get('final_force_maximum')
            if fmax is not None:
                results.final_force_maximum = fmax * force_unit

        sec_workflow = GeometryOptimization()
        sec_workflow.method = model
        sec_workflow.results = results
        self.archive.workflow2 = sec_workflow

    def parse_workflow(self, simulation: Simulation) -> None:
        if not simulation.outputs:
            return
        min_stats = self._log_parser.get('minimization_stats')
        if min_stats is not None:
            self._parse_geometry_optimization_workflow()
        else:
            self.parse_md_workflow({})
            if (
                self._md_method is not None
                and self.archive is not None
                and self.archive.workflow2 is not None
            ):
                self.archive.workflow2.method = self._md_method

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
            _atom_style_columns = {
                'atomic': 'id type x y z',
                'charge': 'id type q x y z',
                'full': 'id mol-id type q x y z',
                'bond': 'id mol-id type x y z',
                'molecular': 'id mol-id type x y z',
                'angle': 'id mol-id type x y z',
                'sphere': 'id type diameter density x y z',
                'dipole': 'id type q x y z mux muy muz',
            }
            atom_style_raw = self._log_parser.get('atom_style')
            atom_style_key = (
                atom_style_raw[0]
                if isinstance(atom_style_raw, list)
                else atom_style_raw
            )
            atom_style_cols = _atom_style_columns.get(atom_style_key)
            kwargs = {'atom_style': atom_style_cols} if atom_style_cols else {}
            traj_parser = MDAnalysisParser(
                topology_format='DATA', format=file_type.upper(), **kwargs
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
            return _create_formatted_parser(traj_file, file_type, data_file)

        # TODO: 'atom' keyword is a LB edit, test
        elif file_type == 'custom' or file_type == 'atom' and data_file:
            return _create_custom_parser(traj_file, index, data_file, file_type)

        elif file_type == 'cfg':
            self.logger.warning(
                'AtomEye CFG dump format is not currently supported; '
                'system section will be empty.'
            )
            traj_parser = TrajParser()
            traj_parser.mainfile = traj_file
            return traj_parser

        else:
            self.logger.warning('File type not recognized.', file_type=traj_file)
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
        # self.parse_input(self.archive.data)

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
        if not parsers:
            self.logger.warning(
                'No trajectory parsers could be created. '
                'Check that trajectory/dump files are present alongside the log file.'
            )
            return
        if self.traj_parsers[0] is None:
            return

        # Parse system, method, parameters, thermodynamic data, etc.
        self._parse_content_sections()

        if not self.archive.data.model_system:
            expected = []
            if self._data_parser.mainfile:
                expected.append(os.path.basename(self._data_parser.mainfile))
            for p in self.traj_parsers._parsers:
                mf = getattr(p, 'mainfile', None)
                if mf:
                    expected.append(os.path.basename(mf))
            self.logger.warning(
                'No system representation could be created. '
                'Upload the referenced structural files to enable '
                'materials recognition.',
                expected_files=expected,
            )

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
