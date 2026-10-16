#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
from ase.data import atomic_numbers
from ase.symbols import symbols2numbers
from nomad.units import ureg
from nomad_simulations.schema_packages.atoms_state import (
    AtomsState,
    CGBeadState,
    ParticleState,
)


def particle_parameters_from_labels(
    labels: Iterable[str | None] | None,
) -> list[dict[str, Any]]:
    """Create particle identity parameters from otherwise unqualified labels.

    Label-only formats do not provide enough information to classify individual
    particles independently. Preserve their established collection-wide policy:
    labels describe atoms only when every label is a valid element symbol;
    otherwise they all describe coarse-grained beads.
    """
    label_list = [] if labels is None else list(labels)
    if not label_list:
        return []

    is_atomic = False
    if all(label not in (None, '') for label in label_list):
        try:
            symbols2numbers(label_list)
            is_atomic = True
        except (KeyError, TypeError, ValueError):
            pass

    return [
        {
            'label': label,
            **({'chemical_symbol': label} if is_atomic else {}),
        }
        for label in label_list
    ]


def _valid_mass(mass: Any) -> Any | None:
    """Return a positive finite mass with units, or ``None`` when invalid."""
    if mass is None:
        return None
    try:
        magnitude = mass.magnitude if hasattr(mass, 'magnitude') else mass
        if not np.isfinite(magnitude) or magnitude <= 0:
            return None
        return mass if hasattr(mass, 'units') else float(mass) * ureg.amu
    except (TypeError, ValueError):
        return None


def particle_state_payloads_from_parameters(
    parameters: Iterable[Mapping[str, Any]] | None,
) -> list[dict[str, Any]]:
    """Build polymorphic particle-state payloads from parser-provided evidence.

    ``chemical_symbol`` must represent an element identity established by the
    parser, rather than a guess from ``label``. Invalid or absent symbols produce
    :class:`CGBeadState` payloads. Numeric masses are interpreted as atomic mass
    units; unit-bearing masses retain their units. Other source properties, such
    as force-field partial charges, are intentionally not copied into intrinsic
    particle-state fields.
    """
    payloads: list[dict[str, Any]] = []
    if parameters is None:
        return payloads

    for parameter in parameters:
        label = parameter.get('label')
        label = str(label) if label is not None else None
        chemical_symbol = parameter.get('chemical_symbol')
        chemical_symbol = str(chemical_symbol) if chemical_symbol is not None else None

        if chemical_symbol in atomic_numbers:
            payload: dict[str, Any] = {
                'chemical_symbol': chemical_symbol,
                'label': label,
            }
        else:
            payload = {
                'm_def': CGBeadState.m_def.qualified_name(),
                'label': label,
            }
            if label:
                payload['bead_symbol'] = str(label)

        mass = _valid_mass(parameter.get('mass'))
        if mass is not None:
            payload['mass'] = mass
        payloads.append(payload)

    return payloads


def particle_states_from_parameters(
    parameters: Iterable[Mapping[str, Any]] | None,
) -> list[ParticleState]:
    """Instantiate particle states from parser-provided identity parameters."""
    particle_states: list[ParticleState] = []
    for payload in particle_state_payloads_from_parameters(parameters):
        common = {'label': payload.get('label')}
        if payload.get('mass') is not None:
            common['mass'] = payload['mass']
        if 'chemical_symbol' in payload:
            particle_states.append(
                AtomsState(
                    chemical_symbol=payload.get('chemical_symbol'),
                    **common,
                )
            )
        else:
            particle_states.append(
                CGBeadState(
                    bead_symbol=payload.get('bead_symbol'),
                    **common,
                )
            )
    return particle_states
