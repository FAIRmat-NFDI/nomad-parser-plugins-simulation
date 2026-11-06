"""Extract MD method settings from the ordered commands echoed in LAMMPS logs."""

from collections.abc import Iterable
from typing import Any

import numpy as np
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    BarostatParameters,
    MolecularDynamicsMethod,
    ThermostatParameters,
)


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
    style = fix[2].lower()
    thermostat_type = {
        'nvt': 'nose_hoover',
        'npt': 'nose_hoover',
        'langevin': 'langevin_schneider',
        'temp/berendsen': 'berendsen',
    }.get(style)
    if thermostat_type is None:
        return
    temperature = (
        _control_values(fix, 'temp')
        if style in ('nvt', 'npt')
        else tuple(_number(value) for value in fix[3:6])
    )
    if len(temperature) != 3 or any(value is None for value in temperature):  # noqa: PLR2004
        return
    start, end, damping = temperature
    thermostat = ThermostatParameters(thermostat_type=thermostat_type)
    if units.get('temperature') is not None:
        if start == end:
            thermostat.temperature_profile = 'constant'
            thermostat.reference_temperature = start * units['temperature']
        else:
            thermostat.temperature_profile = 'linear'
            thermostat.reference_temperature_start = start * units['temperature']
            thermostat.reference_temperature_end = end * units['temperature']
    if units.get('time') is not None:
        thermostat.coupling_constant = damping * units['time']
    method.thermostat_parameters.append(thermostat)


def _add_pressure_control(
    method: MolecularDynamicsMethod, fix: list[str], units: dict[str, Any]
) -> None:
    if fix[2].lower() not in ('npt', 'nph'):
        return
    barostat = BarostatParameters(barostat_type='nose_hoover')
    pressure = _control_values(fix, 'iso')
    if pressure is not None:
        start, end, damping = pressure
        if units.get('pressure') is not None:
            if start == end:
                barostat.pressure_profile = 'constant'
                barostat.reference_pressure = np.eye(3) * start * units['pressure']
            else:
                barostat.pressure_profile = 'linear'
                barostat.reference_pressure_start = (
                    np.eye(3) * start * units['pressure']
                )
                barostat.reference_pressure_end = np.eye(3) * end * units['pressure']
        if units.get('time') is not None:
            barostat.coupling_constant = np.eye(3) * damping * units['time']
    method.barostat_parameters.append(barostat)


def build_md_method(
    settings: dict[str, Any], units: dict[str, Any]
) -> MolecularDynamicsMethod:
    """Map one selected LAMMPS run to the shared MD method schema."""
    method = MolecularDynamicsMethod()
    fixes = settings['fixes']
    styles = [
        str(fix[2]).lower()
        for fix in fixes
        if len(fix) >= 3 and fix[1] == 'all'  # noqa: PLR2004
    ]
    integrators = [style for style in styles if style in ('nve', 'nvt', 'npt', 'nph')]
    has_langevin = 'langevin' in styles
    has_berendsen = 'temp/berendsen' in styles
    if len(integrators) == 1:
        ensemble = integrators[0].upper()
        if ensemble == 'NVE' and (has_langevin or has_berendsen):
            ensemble = 'NVT'
        method.thermodynamic_ensemble = ensemble
        if settings['run_style'] is None or settings['run_style'] == 'verlet':
            method.integrator_type = 'velocity_verlet'
        elif settings['run_style'] == 'respa':
            method.integrator_type = 'rRESPA_multitimescale'

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
        if len(fix) < 4 or fix[1] != 'all':  # noqa: PLR2004
            continue
        _add_temperature_control(method, fix, units)
        _add_pressure_control(method, fix, units)
    return method
