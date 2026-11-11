"""Source-backed LAMMPS force-field settings and simple potentials."""

from collections import defaultdict
from typing import Any

import numpy as np
from nomad.units import ureg
from nomad_simulations.schema_packages.force_field import (
    ForceCalculations,
    ForceField,
    HarmonicAngle,
    HarmonicBond,
    ParameterEntry,
    ParticleParameters,
    ParticleParametersContainer,
    Potential,
)

HARMONIC_COEFFICIENT_COUNT = 2
MIN_PAIR_COEFFICIENT_COUNT = 4
MAX_PAIR_VALUES = 4
TABULAR_DIMENSIONS = 2
MASS_ROW_LENGTH = 2


def _last_command(reader: Any, name: str) -> list[Any] | None:
    commands = reader.get(name)
    if not commands:
        return None
    command = commands[-1]
    if isinstance(command, (list, tuple, np.ndarray)):
        return list(command)
    return [command]


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _coefficients(log: Any, data: Any, kind: str) -> dict[int, list[float]]:
    """Return explicit numeric type coefficients, with log commands overriding data."""
    coefficients: dict[int, list[float]] = {}
    if data.mainfile is not None:
        for _, rows in data.get(f'{kind.title()} Coeffs') or []:
            if rows is None:
                continue
            for row in rows:
                values = [_number(value) for value in row]
                if any(value is None for value in values) or not values:
                    continue
                type_id = values[0]
                if type_id.is_integer() and type_id > 0:
                    coefficients[int(type_id)] = values[1:]
    for row in log.get(f'{kind}_coeff') or []:
        values = [_number(value) for value in row]
        if any(value is None for value in values) or not values:
            continue
        type_id = values[0]
        if type_id.is_integer() and type_id > 0:
            coefficients[int(type_id)] = values[1:]
    return coefficients


def _typed_rows(  # noqa: PLR0911
    data: Any, kind: str, atom_ids: Any
) -> dict[int, list[list[int]]]:
    """Map data-file atom IDs to zero-based trajectory positions when exact."""
    if data.mainfile is None or atom_ids is None:
        return {}
    sections = data.get(f'{kind.title()}s') or []
    if not sections or sections[0][1] is None:
        return {}
    ids = np.asarray(atom_ids)
    if ids.ndim != 1 or not np.isfinite(ids).all():
        return {}
    if not np.array_equal(ids, ids.astype(int)):
        return {}
    id_to_index = {int(atom_id): index for index, atom_id in enumerate(ids)}
    if len(id_to_index) != len(ids):
        return {}

    width = {'bond': 2, 'angle': 3, 'dihedral': 4}[kind]
    result: dict[int, list[list[int]]] = defaultdict(list)
    for _, rows in sections:
        if rows is None:
            continue
        for row in rows:
            if len(row) != width + 2 or not np.isfinite(row).all():
                return {}
            if not np.array_equal(row, row.astype(int)):
                return {}
            type_id = int(row[1])
            try:
                result[type_id].append(
                    [id_to_index[int(atom_id)] for atom_id in row[2:]]
                )
            except KeyError:
                return {}
    return result


def _force_calculations(  # noqa: PLR0912
    log: Any, data: Any, units: dict[str, Any]
) -> ForceCalculations | None:
    settings = ForceCalculations()
    style = _last_command(log, 'pair_style')
    if style:
        pair_style = str(style[0]).lower()
        args = [_number(value) for value in style[1:]]
        if pair_style == 'lj/cut' and len(args) == 1 and args[0] is not None:
            settings.vdw_cutoff = args[0] * units['distance']
        elif pair_style in {'lj/cut/coul/cut', 'lj/cut/coul/long'}:
            if len(args) in (1, 2) and all(value is not None for value in args):
                settings.vdw_cutoff = args[0] * units['distance']
                settings.coulomb_cutoff = args[-1] * units['distance']
        elif pair_style == 'coul/cut' and len(args) == 1 and args[0] is not None:
            settings.coulomb_cutoff = args[0] * units['distance']

    kspace = _last_command(log, 'kspace_style')
    if kspace:
        kspace_style = str(kspace[0]).lower()
        settings.coulomb_type = {
            'ewald': 'ewald',
            'pppm': 'particle_particle_particle_mesh',
            'msm': 'multilevel_summation',
        }.get(kspace_style)
    elif settings.coulomb_cutoff is not None:
        settings.coulomb_type = 'cutoff'

    neighbor = _last_command(log, 'neighbor')
    cutoffs = [
        cutoff
        for cutoff in (settings.vdw_cutoff, settings.coulomb_cutoff)
        if cutoff is not None
    ]
    if style and str(style[0]).lower() in {
        'lj/cut',
        'lj/cut/coul/cut',
        'lj/cut/coul/long',
    }:
        pair_style = str(style[0]).lower()
        allowed_lengths = (2, 3, 4) if pair_style == 'lj/cut/coul/cut' else (2, 3)
        for values in _pair_coefficients(log, data).values():
            if len(values) in allowed_lengths:
                cutoffs.extend(value * units['distance'] for value in values[2:])
    if neighbor and cutoffs:
        skin = _number(neighbor[0])
        if skin is not None:
            settings.neighbor_update_cutoff = max(cutoffs) + skin * units['distance']

    neighbor_options = {}
    for command in log.get('neigh_modify') or []:
        tokens = [str(value).lower() for value in command]
        for index, token in enumerate(tokens[:-1]):
            if token in {'every', 'check', 'once'}:
                neighbor_options[token] = tokens[index + 1]
    # With check=yes, ``every`` is only a consideration interval, not an
    # actual rebuild frequency; once=yes suppresses subsequent rebuilding.
    if neighbor_options.get('check') == 'no' and neighbor_options.get('once') != 'yes':
        every = _number(neighbor_options.get('every'))
        if every is not None and every.is_integer() and every > 0:
            settings.neighbor_update_frequency = int(every)

    return settings if settings.m_to_dict() else None


def _harmonic_potentials(
    log: Any, data: Any, atom_ids: Any, kind: str, units: dict[str, Any]
) -> list[Potential]:
    style = _last_command(log, f'{kind}_style')
    if not style or str(style[0]).lower() != 'harmonic':
        return []
    coefficients = _coefficients(log, data, kind)
    rows_by_type = _typed_rows(data, kind, atom_ids)
    result = []
    for type_id, values in sorted(coefficients.items()):
        if len(values) != HARMONIC_COEFFICIENT_COUNT:
            continue
        constant, equilibrium = values
        if kind == 'bond':
            potential = HarmonicBond(
                type='bond',
                functional_form='harmonic',
                force_constant=2 * constant * units['energy'] / units['distance'] ** 2,
                equilibrium_value=equilibrium * units['distance'],
            )
        else:
            potential = HarmonicAngle(
                type='angle',
                functional_form='harmonic',
                force_constant=2 * constant * units['energy'] / ureg.radian**2,
                equilibrium_value=equilibrium * ureg.degree,
            )
        potential.parameters.append(
            ParameterEntry(name='lammps_type', value=str(type_id))
        )
        if rows := rows_by_type.get(type_id):
            potential.particle_indices = np.asarray(rows, dtype=np.int32)
        result.append(potential)
    return result


def _pair_coefficients(log: Any, data: Any) -> dict[tuple[int, int], list[float]]:
    """Collect explicit LJ coefficients; command values override data values."""
    coefficients: dict[tuple[int, int], list[float]] = {}
    if data.mainfile is not None:
        for section, n_type_columns in (('Pair Coeffs', 1), ('PairIJ Coeffs', 2)):
            for _, rows in data.get(section) or []:
                if rows is None:
                    continue
                for row in rows:
                    values = [_number(value) for value in row]
                    if any(value is None for value in values):
                        continue
                    if len(values) < n_type_columns + 2:
                        continue
                    types = values[:n_type_columns]
                    if not all(value.is_integer() and value > 0 for value in types):
                        continue
                    pair = (int(types[0]), int(types[-1]))
                    coefficients[tuple(sorted(pair))] = values[n_type_columns:]
    for row in log.get('pair_coeff') or []:
        values = [_number(value) for value in row]
        if len(values) < MIN_PAIR_COEFFICIENT_COUNT or any(
            value is None for value in values
        ):
            continue
        type_i, type_j = values[:2]
        if not (type_i.is_integer() and type_j.is_integer()):
            continue
        coefficients[tuple(sorted((int(type_i), int(type_j))))] = values[2:]
    return coefficients


def _pair_potentials(log: Any, data: Any, units: dict[str, Any]) -> list[Potential]:
    style = _last_command(log, 'pair_style')
    if not style or str(style[0]).lower() not in {
        'lj/cut',
        'lj/cut/coul/cut',
        'lj/cut/coul/long',
    }:
        return []
    result = []
    for (type_i, type_j), values in sorted(_pair_coefficients(log, data).items()):
        if len(values) not in (2, 3, 4):
            continue
        if (
            len(values) == MAX_PAIR_VALUES
            and str(style[0]).lower() != 'lj/cut/coul/cut'
        ):
            continue
        epsilon, sigma = values[:2]
        potential = Potential(
            name=f'{style[0]} types {type_i}-{type_j}',
            type='nonbonded',
            functional_form='Lennard-Jones 12-6',
        )
        for name, value, unit in (
            ('type_i', str(type_i), None),
            ('type_j', str(type_j), None),
            (
                'epsilon',
                str((epsilon * units['energy']).to('joule').magnitude),
                'joule',
            ),
            ('sigma', str((sigma * units['distance']).to('meter').magnitude), 'meter'),
        ):
            potential.parameters.append(
                ParameterEntry(name=name, value=value, unit=unit)
            )
        for index, cutoff in enumerate(values[2:]):
            name = 'vdw_cutoff' if index == 0 else 'coulomb_cutoff'
            potential.parameters.append(
                ParameterEntry(
                    name=name,
                    value=str((cutoff * units['distance']).to('meter').magnitude),
                    unit='meter',
                )
            )
        # Pair contacts change during a trajectory; coefficient declarations
        # are not a count of actual nonbonded interactions.
        result.append(potential)
    return result


def _opls_potentials(
    log: Any, data: Any, atom_ids: Any, units: dict[str, Any]
) -> list[Potential]:
    style = _last_command(log, 'dihedral_style')
    if not style or str(style[0]).lower() != 'opls':
        return []
    coefficients = _coefficients(log, data, 'dihedral')
    rows_by_type = _typed_rows(data, 'dihedral', atom_ids)
    result = []
    for type_id, values in sorted(coefficients.items()):
        if len(values) != MAX_PAIR_VALUES:
            continue
        potential = Potential(
            name=f'OPLS dihedral type {type_id}',
            type='dihedral',
            functional_form='OPLS',
        )
        potential.parameters.append(
            ParameterEntry(name='lammps_type', value=str(type_id))
        )
        for index, value in enumerate(values, start=1):
            potential.parameters.append(
                ParameterEntry(
                    name=f'K{index}',
                    value=str((value * units['energy']).to('joule').magnitude),
                    unit='joule',
                )
            )
        if rows := rows_by_type.get(type_id):
            potential.particle_indices = np.asarray(rows, dtype=np.int32)
        result.append(potential)
    return result


def build_force_field(
    log: Any,
    data: Any,
    atom_ids: Any,
    units: dict[str, Any],
    logger: Any = None,
) -> ForceField | None:
    """Build only force-field terms supported by explicit, unambiguous input."""
    supported_styles = {
        'pair': {'lj/cut', 'lj/cut/coul/cut', 'lj/cut/coul/long', 'coul/cut'},
        'bond': {'harmonic'},
        'angle': {'harmonic'},
        'dihedral': {'opls'},
        'improper': set(),
    }
    if logger is not None:
        for kind, supported in supported_styles.items():
            style = _last_command(log, f'{kind}_style')
            if style and str(style[0]).lower() not in supported | {'none'}:
                logger.warning(
                    f'LAMMPS {kind}_style {style[0]} is not mapped to '
                    'force-field contributions.'
                )
    settings = _force_calculations(log, data, units)
    contributions = _pair_potentials(log, data, units)
    if logger is not None and log.get('pair_coeff') and not contributions:
        logger.warning(
            'LAMMPS pair coefficients could not be mapped as explicit numeric '
            'Lennard-Jones type pairs.'
        )
    for kind in ('bond', 'angle'):
        terms = _harmonic_potentials(log, data, atom_ids, kind, units)
        if logger is not None and log.get(f'{kind}_coeff') and not terms:
            logger.warning(
                f'LAMMPS {kind} coefficients could not be mapped as explicit '
                'numeric harmonic types.'
            )
        contributions.extend(terms)
    dihedrals = _opls_potentials(log, data, atom_ids, units)
    if logger is not None and log.get('dihedral_coeff') and not dihedrals:
        logger.warning(
            'LAMMPS dihedral coefficients could not be mapped as explicit '
            'numeric OPLS types.'
        )
    contributions.extend(dihedrals)
    if settings is None and not contributions:
        return None
    force_field = ForceField()
    if settings is not None:
        force_field.numerical_settings.append(settings)
    force_field.contributions.extend(contributions)
    return force_field


def attach_particle_parameters(  # noqa: PLR0913, PLR0917, PLR0912
    force_field: ForceField,
    log: Any,
    data: Any,
    atom_ids: Any,
    particle_states: Any,
    units: dict[str, Any],
) -> None:
    """Attach data-file masses and partial charges to exact particle refs."""
    if data.mainfile is None or atom_ids is None or not particle_states:
        return
    style_command = _last_command(log, 'atom_style')
    style = str(style_command[0]).lower() if style_command else None
    columns = {
        'atomic': (1, None, 5),
        'charge': (1, 2, 6),
        'molecular': (2, None, 6),
        'full': (2, 3, 7),
    }
    if style not in columns:
        return
    sections = data.get('Atoms') or []
    if not sections or sections[0][1] is None:
        return
    atoms = sections[0][1]
    type_column, charge_column, width = columns[style]
    ids = np.asarray(atom_ids)
    if (
        atoms.ndim != TABULAR_DIMENSIONS
        or atoms.shape[1] not in (width, width + 3)
        or len(atoms) != len(ids)
        or len(ids) != len(particle_states)
        or not np.isfinite(atoms).all()
        or not np.isfinite(ids).all()
        or not np.array_equal(atoms[:, 0], atoms[:, 0].astype(int))
        or not np.array_equal(atoms[:, type_column], atoms[:, type_column].astype(int))
        or not np.array_equal(ids, ids.astype(int))
    ):
        return
    id_to_index = {int(atom_id): index for index, atom_id in enumerate(ids)}
    if len(id_to_index) != len(ids) or set(atoms[:, 0].astype(int)) != set(id_to_index):
        return

    masses = {}
    for _, rows in data.get('Masses') or []:
        if rows is None:
            continue
        for row in rows:
            if len(row) == MASS_ROW_LENGTH and row[0].is_integer() and row[1] > 0:
                masses[int(row[0])] = float(row[1])

    grouped: dict[tuple[int, float | None], list[Any]] = defaultdict(list)
    for row in atoms:
        type_id = int(row[type_column])
        charge = float(row[charge_column]) if charge_column is not None else None
        grouped[(type_id, charge)].append(particle_states[id_to_index[int(row[0])]])

    container = ParticleParametersContainer()
    for (type_id, charge), states in sorted(grouped.items()):
        if type_id not in masses and charge is None:
            continue
        parameters = ParticleParameters(
            particle_type=str(type_id),
            species_scope=states,
        )
        if type_id in masses:
            parameters.effective_mass = masses[type_id] * units['mass']
        if charge is not None:
            parameters.partial_charge = charge * units['charge']
        container.particle_parameters.append(parameters)
    if container.particle_parameters:
        force_field.numerical_settings.append(container)
