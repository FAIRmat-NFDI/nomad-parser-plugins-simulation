"""Extract MD workflow-method settings from ordered LAMMPS log commands."""

from collections.abc import Iterable
from typing import Any

import numpy as np
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    BarostatParameters,
    MolecularDynamicsMethod,
    ThermostatParameters,
)

N_CONTROL_VALUES = 3


def _number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _frequency(value: str) -> int | None:
    number = _number(value)
    if number is not None and number > 0 and number.is_integer():
        return int(number)
    return None


def _fix_style(fix: list[str]) -> str:
    """Ignore LAMMPS accelerator suffixes, which do not change fix semantics."""
    style = fix[2].lower()
    base, separator, suffix = style.rpartition('/')
    accelerated = separator and suffix in ('gpu', 'intel', 'kk', 'omp', 'opt')
    return base if accelerated else style


def selected_run(lines: Iterable[str]) -> dict[str, Any] | None:  # noqa: PLR0912
    """Return settings active at the last run or minimize command.

    The log normally echoes commands in execution order. An auxiliary log is
    read after the main log, so settings declared before ``log`` remain active.
    Variable expressions are skipped until an expanded command is echoed.
    """
    fixes: dict[str, list[str]] = {}
    dumps: dict[str, list[str | int | None]] = {}
    timestep = thermo = run_style = None
    current_step: int | None = None
    selected = None

    for line in lines:
        command = line.split('#', 1)[0].split()
        if not command or any('$' in token for token in command):
            continue
        name, *args = command
        name = name.lower()
        if name == 'fix' and len(args) >= 3:  # noqa: PLR2004
            fixes[args[0]] = args
        elif name == 'unfix' and args:
            fixes.pop(args[0], None)
        elif name == 'dump' and len(args) >= 5:  # noqa: PLR2004
            dumps[args[0]] = [*args, _frequency(args[3])]
        elif name == 'undump' and args:
            dumps.pop(args[0], None)
        elif name == 'dump_modify' and len(args) >= 3 and args[0] in dumps:  # noqa: PLR2004
            if 'every' in args[1:]:
                i = args.index('every')
                dumps[args[0]][-1] = (
                    _frequency(args[i + 1]) if i + 1 < len(args) else None
                )
            elif 'every/time' in args[1:]:
                dumps[args[0]][-1] = None
        elif name == 'timestep' and args:
            timestep = _number(args[0])
        elif name == 'thermo' and args:
            thermo = _frequency(args[0])
        elif name == 'run_style' and args:
            run_style = args[0].lower()
        elif name == 'reset_timestep' and args:
            value = _number(args[0])
            current_step = (
                int(value) if value is not None and value.is_integer() else None
            )
        elif name == 'run' and args:
            value = _number(args[0])
            steps = (
                int(value)
                if value is not None and value.is_integer() and value >= 0
                else None
            )
            if 'upto' in args[1:]:
                steps = (
                    max(0, steps - current_step)
                    if steps is not None and current_step is not None
                    else None
                )
            selected = dict(
                kind='md',
                fixes=list(fixes.values()),
                dumps=[list(d) for d in dumps.values()],
                timestep=timestep,
                thermo=thermo,
                run_style=run_style,
                n_steps=steps,
            )
            current_step = (
                current_step + steps
                if current_step is not None and steps is not None
                else None
            )
        elif name in ('minimize', 'minimize/kk'):
            selected = {'kind': 'minimize'}
            # Minimization can advance the step counter by an unknown amount.
            current_step = None
    return selected


def _one_frequency(dumps: list[list[str | int | None]], kind: str) -> int | None:
    frequencies = set()
    for dump in dumps:
        style = str(dump[2]).lower()
        fields = set(str(field).lower() for field in dump[5:-1])
        if kind == 'coordinate':
            has_kind = style in ('atom', 'xyz', 'dcd') or (
                style == 'custom'
                and all(
                    any(f in fields for f in variants)
                    for variants in (
                        ('x', 'xu', 'xs', 'xsu'),
                        ('y', 'yu', 'ys', 'ysu'),
                        ('z', 'zu', 'zs', 'zsu'),
                    )
                )
            )
        else:
            prefix = 'v' if kind == 'velocity' else 'f'
            has_kind = style == 'custom' and all(
                f'{prefix}{axis}' in fields for axis in 'xyz'
            )
        if has_kind:
            if not isinstance(dump[-1], int):
                return None
            frequencies.add(dump[-1])
    return frequencies.pop() if len(frequencies) == 1 else None


def _control_values(fix: list[str], keyword: str) -> tuple[float, float, float] | None:
    if keyword not in fix[3:]:
        return None
    i = fix.index(keyword, 3)
    if i + 3 >= len(fix):
        return None
    values = tuple(_number(value) for value in fix[i + 1 : i + 4])
    return values if all(value is not None for value in values) else None


def _add_temperature_control(
    method: MolecularDynamicsMethod, fix: list[str], units: dict[str, Any]
) -> None:
    style = _fix_style(fix)
    thermostat_type = {
        'nvt': 'nose_hoover',
        'npt': 'nose_hoover',
        'langevin': 'langevin_schneider',
        'temp/berendsen': 'berendsen',
        'temp/csvr': 'velocity_rescaling',
        'temp/csld': 'velocity_rescaling_langevin',
        'brownian': 'brownian',
        'brownian/sphere': 'brownian',
        'brownian/asphere': 'brownian',
    }.get(style)
    if thermostat_type is None:
        return
    thermostat = ThermostatParameters(thermostat_type=thermostat_type)
    if style.startswith('brownian'):
        # Brownian fixes take one temperature and a seed; gamma_t is a
        # friction coefficient, not a temperature-coupling time constant.
        temperature = _number(fix[3]) if len(fix) > 3 else None  # noqa: PLR2004
        if temperature is not None and units.get('temperature') is not None:
            thermostat.temperature_profile = 'constant'
            thermostat.reference_temperature = temperature * units['temperature']
    else:
        temperature = (
            _control_values(fix, 'temp')
            if style in ('nvt', 'npt')
            else tuple(_number(value) for value in fix[3:6])
        )
        if (
            temperature is not None
            and len(temperature) == N_CONTROL_VALUES
            and all(value is not None for value in temperature)
        ):
            start, end, damping = temperature
            if units.get('temperature') is not None:
                if start == end:
                    thermostat.temperature_profile = 'constant'
                    thermostat.reference_temperature = start * units['temperature']
                else:
                    thermostat.temperature_profile = 'linear'
                    thermostat.reference_temperature_start = (
                        start * units['temperature']
                    )
                    thermostat.reference_temperature_end = end * units['temperature']
            if units.get('time') is not None:
                thermostat.coupling_constant = damping * units['time']
    method.thermostat_parameters.append(thermostat)


def _pressure_components(
    fix: list[str], style: str
) -> tuple[dict[str, tuple[float, float, float]], str | None]:
    """Expand LAMMPS pressure shortcuts and retain explicit tensor components."""
    components: dict[str, tuple[float, float, float]] = {}
    coupling_type = None
    for keyword in ('iso', 'aniso', 'tri'):
        if keyword == 'tri' and style == 'press/berendsen':
            continue
        values = _control_values(fix, keyword)
        if values is not None:
            components.update({axis: values for axis in 'xyz'})
            if keyword == 'tri':
                components.update(
                    {axis: (0.0, 0.0, values[2]) for axis in ('xy', 'xz', 'yz')}
                )
            coupling_type = 'isotropic' if keyword == 'iso' else 'anisotropic'

    axes = (
        ('x', 'y', 'z')
        if style == 'press/berendsen'
        else ('x', 'y', 'z', 'xy', 'xz', 'yz')
    )
    for axis in axes:
        values = _control_values(fix, axis)
        if values is not None:
            components[axis] = values

    if 'couple' in fix[3:]:
        i = fix.index('couple', 3)
        if i + 1 < len(fix):
            coupling_type = {
                'xyz': 'isotropic',
                'xy': 'semi_isotropic',
                'xz': 'semi_isotropic',
                'yz': 'semi_isotropic',
                'none': 'anisotropic',
            }.get(fix[i + 1].lower())
    elif components and coupling_type is None:
        coupling_type = 'anisotropic'
    return components, coupling_type


def _pressure_matrices(
    components: dict[str, tuple[float, float, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Represent controlled pressure directions in symmetric 3×3 arrays."""

    matrix_positions = {
        'x': ((0, 0),),
        'y': ((1, 1),),
        'z': ((2, 2),),
        'xy': ((0, 1), (1, 0)),
        'xz': ((0, 2), (2, 0)),
        'yz': ((1, 2), (2, 1)),
    }
    starts, ends, damping = (np.zeros((3, 3)) for _ in range(3))
    for axis, values in components.items():
        for row, column in matrix_positions[axis]:
            starts[row, column], ends[row, column], damping[row, column] = values
    return starts, ends, damping


def _add_pressure_control(
    method: MolecularDynamicsMethod, fix: list[str], units: dict[str, Any]
) -> None:
    style = _fix_style(fix)
    if style not in ('npt', 'nph', 'press/berendsen'):
        return
    barostat = BarostatParameters(
        barostat_type='berendsen' if style == 'press/berendsen' else 'nose_hoover'
    )
    components, coupling_type = _pressure_components(fix, style)
    barostat.coupling_type = coupling_type
    starts, ends, damping = _pressure_matrices(components)

    # An absent target axis is not a target of zero pressure. Keep the target
    # matrix unset unless all three diagonal components are known.
    pressure_unit = units.get('pressure')
    if pressure_unit is not None and all(axis in components for axis in 'xyz'):
        if np.array_equal(starts, ends):
            barostat.pressure_profile = 'constant'
            barostat.reference_pressure = starts * pressure_unit
        else:
            barostat.pressure_profile = 'linear'
            barostat.reference_pressure_start = starts * pressure_unit
            barostat.reference_pressure_end = ends * pressure_unit
    if units.get('time') is not None and components:
        barostat.coupling_constant = damping * units['time']

    if style == 'press/berendsen' and 'modulus' in fix[3:]:
        i = fix.index('modulus', 3)
        modulus = _number(fix[i + 1]) if i + 1 < len(fix) else None
        if modulus is not None and modulus > 0 and pressure_unit is not None:
            barostat.compressibility = np.eye(3) / (modulus * pressure_unit)
    method.barostat_parameters.append(barostat)


def _set_integrator_and_ensemble(
    method: MolecularDynamicsMethod, styles: list[str], run_style: str | None
) -> None:
    """Describe a unique all-particle integrator and its active controls."""
    integrators = [
        style
        for style in styles
        if style
        in (
            'nve',
            'nvt',
            'npt',
            'nph',
            'brownian',
            'brownian/sphere',
            'brownian/asphere',
        )
    ]
    if len(integrators) != 1:
        return

    integrator = integrators[0]
    temperature_styles = (
        'nvt',
        'npt',
        'langevin',
        'temp/berendsen',
        'temp/csvr',
        'temp/csld',
    )
    has_temperature = any(style in temperature_styles for style in styles)
    has_pressure = any(style in ('npt', 'nph', 'press/berendsen') for style in styles)
    if integrator.startswith('brownian'):
        method.thermodynamic_ensemble = 'NPT' if has_pressure else 'NVT'
        method.integrator_type = 'brownian'
        return

    ensemble = integrator.upper()
    if integrator == 'nve':
        if has_temperature and has_pressure:
            ensemble = 'NPT'
        elif has_temperature:
            ensemble = 'NVT'
        elif has_pressure:
            ensemble = 'NPH'
    elif (integrator == 'nvt' and has_pressure) or (
        integrator == 'nph' and has_temperature
    ):
        ensemble = 'NPT'
    method.thermodynamic_ensemble = ensemble
    if run_style is None or run_style == 'verlet':
        method.integrator_type = 'velocity_verlet'
    elif run_style == 'respa':
        method.integrator_type = 'rRESPA_multitimescale'


def build_md_method(
    settings: dict[str, Any], units: dict[str, Any]
) -> MolecularDynamicsMethod:
    """Map one selected LAMMPS run to the shared MD method schema."""
    method = MolecularDynamicsMethod()
    fixes = settings['fixes']
    styles = [
        _fix_style(fix)
        for fix in fixes
        if len(fix) >= 3 and fix[1].lower() == 'all'  # noqa: PLR2004
    ]
    _set_integrator_and_ensemble(method, styles, settings['run_style'])

    time_unit = units.get('time')
    if settings['timestep'] is not None and time_unit is not None:
        method.integration_timestep = settings['timestep'] * time_unit
    method.n_steps = settings['n_steps']
    method.thermodynamics_save_frequency = settings['thermo']
    for kind in ('coordinate', 'velocity', 'force'):
        frequency = _one_frequency(settings['dumps'], kind)
        if frequency is not None:
            setattr(method, f'{kind}_save_frequency', frequency)

    for fix in fixes:
        if len(fix) < 3 or fix[1].lower() != 'all':  # noqa: PLR2004
            continue
        _add_temperature_control(method, fix, units)
        _add_pressure_control(method, fix, units)
    return method
