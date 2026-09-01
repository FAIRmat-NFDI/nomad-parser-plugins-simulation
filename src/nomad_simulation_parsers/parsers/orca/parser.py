import re
from importlib import reload
from typing import TYPE_CHECKING, Any

import numpy as np
from nomad.datamodel.metainfo.workflow import Link, Task
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad_file_parser.mapping_parser import MetainfoParser
from nomad_file_parser.mapping_parser import TextParser as MappingTextParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_method import (
    CC,
    HF,
    OrbitalLocalization,
    PerturbationMethod,
)
from nomad_simulations.schema_packages.workflow.general import SerialWorkflow

from nomad_simulation_parsers.schema_packages import orca

from .text_parser import OutReader

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import EntryArchive
    from structlog.stdlib import BoundLogger

CARTESIAN_COORDINATE_LENGTH = 3
CARTESIAN_COORDINATE_STRIDE = 4
ORBITAL_ENERGY_N_COLUMNS = 4
ORBITAL_ENERGY_N_DIMENSIONS = 2
ORBITAL_OCCUPATION_COLUMN = 1
ORBITAL_ENERGY_EV_COLUMN = 3
ORBITAL_RANGE_BOUNDS = 2
AO_ROW_RE = re.compile(r'\d+[A-Z][a-z]?$')


def str_to_cartesian_coordinates(
    value: list[Any],
) -> tuple[list[str], np.ndarray]:
    cleaned = [
        item.replace('>', '') if isinstance(item, str) else item
        for item in value
        if item != '>'
    ]
    symbols = []
    coordinates = []
    for index in range(0, len(cleaned), CARTESIAN_COORDINATE_STRIDE):
        symbol = cleaned[index]
        coordinate = cleaned[index + 1 : index + CARTESIAN_COORDINATE_STRIDE]
        if isinstance(symbol, str) and len(coordinate) == CARTESIAN_COORDINATE_LENGTH:
            symbols.append(symbol)
            coordinates.append(coordinate)
    return symbols, np.asarray(coordinates, dtype=np.float64) * ureg.angstrom


def str_to_mo_coefficients(value: str | list[str] | None) -> np.ndarray | None:
    if not value:
        return None

    if isinstance(value, list):
        return _token_list_to_mo_coefficients([str(token) for token in value])

    coefficients_by_mo: dict[int, list[float]] = {}
    mo_indices: list[int] = []
    for line in value.splitlines():
        parts = line.split()
        if not parts:
            continue

        if all(part.isdigit() for part in parts):
            mo_indices = [int(part) for part in parts]
            continue

        if (
            not mo_indices
            or not AO_ROW_RE.fullmatch(parts[0])
            or len(parts) < len(mo_indices) + 2
        ):
            continue

        row_values = [float(part) for part in parts[-len(mo_indices) :]]
        for mo_index, coefficient in zip(mo_indices, row_values):
            coefficients_by_mo.setdefault(mo_index, []).append(coefficient)

    return _coefficient_matrix(coefficients_by_mo)


def _token_list_to_mo_coefficients(tokens: list[str]) -> np.ndarray | None:
    coefficients_by_mo: dict[int, list[float]] = {}
    index = 0

    while index < len(tokens):
        # Keep this as a while loop: each MO block has a variable number of
        # columns and rows, so the start of the next block is only known
        # after the current block has been consumed.
        index = _find_next_mo_header(tokens, index)
        mo_indices, index = _read_mo_indices(tokens, index)
        if not mo_indices:
            break

        n_indices = len(mo_indices)
        index += n_indices * 2
        separator_tokens = tokens[index : index + n_indices]
        index += sum(1 for token in separator_tokens if token.startswith('-'))

        index = _read_coefficient_rows(tokens, index, mo_indices, coefficients_by_mo)

    return _coefficient_matrix(coefficients_by_mo)


def _find_next_mo_header(tokens: list[str], start: int) -> int:
    return next(
        (index for index, token in enumerate(tokens[start:], start) if token.isdigit()),
        len(tokens),
    )


def _read_mo_indices(tokens: list[str], start: int) -> tuple[list[int], int]:
    mo_indices = []
    for index, token in enumerate(tokens[start:], start):
        if not token.isdigit():
            return mo_indices, index
        mo_indices.append(int(token))
    return mo_indices, len(tokens)


def _read_coefficient_rows(
    tokens: list[str],
    index: int,
    mo_indices: list[int],
    coefficients_by_mo: dict[int, list[float]],
) -> int:
    n_indices = len(mo_indices)
    row_stride = n_indices + 2
    for row_index in range(index, len(tokens), row_stride):
        if not AO_ROW_RE.fullmatch(tokens[row_index]):
            return row_index

        row_start = row_index + 2
        row_end = row_start + n_indices
        if row_end > len(tokens):
            return row_index

        try:
            row_values = [float(token) for token in tokens[row_start:row_end]]
        except ValueError:
            return row_index

        for mo_index, coefficient in zip(mo_indices, row_values):
            coefficients_by_mo.setdefault(mo_index, []).append(coefficient)

    return len(tokens)


def _coefficient_matrix(
    coefficients_by_mo: dict[int, list[float]],
) -> np.ndarray | None:
    if not coefficients_by_mo:
        return None

    coefficients = [coefficients_by_mo[index] for index in sorted(coefficients_by_mo)]
    if len({len(row) for row in coefficients}) != 1:
        return None

    return np.asarray(coefficients, dtype=np.float64)


class OutParser(MappingTextParser):
    _BASIS_SET_ROLES = {
        'main_basis_set': 'orbital',
        'auxj_basis_set': 'auxiliary_scf',
        'auxjk_basis_set': 'auxiliary_scf',
        'auxc_basis_set': 'auxiliary_post_hf',
    }
    _BASIS_SET_KEYS = {
        **dict.fromkeys(
            ('HF', 'DFT', 'MultireferenceSCF', 'MultireferenceCI'),
            frozenset({'main_basis_set', 'auxj_basis_set', 'auxjk_basis_set'}),
        ),
        **dict.fromkeys(
            ('PerturbationMethod', 'CC', 'MultireferencePT'),
            frozenset({'main_basis_set', 'auxc_basis_set'}),
        ),
    }

    def __init__(self) -> None:
        super().__init__(text_parser=OutReader())
        self._method = None

    def load_file(self) -> OutReader:
        text_parser = super().load_file()
        text_parser.findlazy = False
        return text_parser

    @staticmethod
    def _parser_results(value: Any) -> dict[str, Any]:
        # Nested TextParser quantities expose parsed values through `_results`.
        return value._results if hasattr(value, '_results') else value or {}

    def _navigate(self, source: dict[str, Any], *keys: str) -> dict[str, Any]:
        current = source
        for key in keys:
            current = self._parser_results(current.get(key))
        return current

    @staticmethod
    def _scalar(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(value, (list, tuple)):
            return value[0] if value else None
        return value

    def _get_cartesian_system(self, source: dict[str, Any]) -> tuple[list[str], Any]:
        single_point = self._navigate(source, 'single_point')
        coordinates = single_point.get('cartesian_coordinates', [])
        return str_to_cartesian_coordinates(coordinates) if coordinates else ([], None)

    @staticmethod
    def _build_particle_states(symbols: list[str]) -> list[dict[str, str]]:
        return [{'chemical_symbol': symbol} for symbol in symbols]

    def _get_charge_and_multiplicity(self, source: dict[str, Any]) -> dict[str, int]:
        scf_settings = self._navigate(
            source, 'single_point', 'self_consistent', 'scf_settings'
        )
        result = {}
        total_charge = self._scalar(scf_settings.get('total_charge'))
        if total_charge is not None:
            result['total_charge'] = int(total_charge)
        multiplicity = self._scalar(scf_settings.get('multiplicity'))
        if multiplicity is not None:
            result['total_spin_multiplicity'] = int(multiplicity)
        return result

    def get_atoms(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        symbols, positions = self._get_cartesian_system(src)
        if not symbols:
            return []

        return [
            {
                'is_representative': True,
                'positions': positions,
                'particle_states': self._build_particle_states(symbols),
                **self._get_charge_and_multiplicity(src),
            }
        ]

    @staticmethod
    def _normalize_localization_method(value: Any) -> str | None:
        raw = OutParser._scalar(value)
        if not isinstance(raw, str):
            return None

        raw = raw.upper()
        for key, normalized in {
            'FOSTER-BOYS': 'Foster-Boys',
            'PIPEK-MEZEY': 'Pipek-Mezey',
            'EDMISTON-RUEDENBERG': 'Edmiston-Ruedenberg',
            'AIOPM-NBO': 'AIOPM-NBO',
            'IBO': 'IBO',
            'NBO': 'NBO',
        }.items():
            if key in raw:
                return normalized
        return None

    @staticmethod
    def _hf_reference_form(value: Any) -> str | None:
        reference = OutParser._scalar(value)
        if not isinstance(reference, str):
            return None
        reference = reference.upper()
        return reference if reference in {'RHF', 'UHF', 'ROHF'} else None

    @staticmethod
    def _dft_reference_form(value: Any) -> str | None:
        reference = OutParser._scalar(value)
        if not isinstance(reference, str):
            return None
        return {
            'RHF': 'RKS',
            'UHF': 'UKS',
            'ROHF': 'ROKS',
            'RKS': 'RKS',
            'UKS': 'UKS',
            'ROKS': 'ROKS',
        }.get(reference.upper())

    def _get_scf_settings(self, source: dict[str, Any]) -> dict[str, Any]:
        return self._navigate(source, 'single_point', 'self_consistent', 'scf_settings')

    def _build_xc(self, scf_settings: dict[str, Any]) -> dict[str, Any]:
        xc = {}
        exchange = self._scalar(scf_settings.get('exchange_functional'))
        correlation = self._scalar(
            scf_settings.get('correlation_functional')
            or scf_settings.get('correl_functional')
        )
        functionals = [
            value
            for value in (exchange, correlation)
            if isinstance(value, str) and value
        ]
        functional_key = '+'.join(dict.fromkeys(functionals))
        if functional_key:
            xc['functional_key'] = functional_key
        hf_fraction = self._scalar(scf_settings.get('fraction_hf_exchange'))
        if hf_fraction is not None:
            xc['global_exact_exchange'] = hf_fraction
        return xc

    def get_dft(self, scf_settings: dict[str, Any]) -> dict[str, Any]:
        xc = self._build_xc(scf_settings)
        if not xc:
            return {}

        method = {'xc': xc}
        reference_form = self._dft_reference_form(scf_settings.get('hf_type'))
        if reference_form:
            method['reference_form'] = reference_form
        return method

    def get_dft_methods(
        self, scf_settings: dict[str, Any]
    ) -> list[dict[str, Any]]:
        method = self.get_dft(scf_settings)
        if method:
            self._method = 'DFT'
        return [method] if method else []

    def _build_active_space_data(self, casscf: dict[str, Any]) -> dict[str, Any]:
        return {
            'n_active_electrons': self._scalar(casscf.get('n_active_electrons')),
            'n_active_orbitals': self._scalar(casscf.get('n_active_orbitals')),
            'orbital_space_type': 'CAS',
        }

    def _collect_state_data(self, casscf: dict[str, Any]) -> dict[str, Any]:
        state_multiplicities = []
        n_roots_per_multiplicity = []
        state_weights = []
        for block in casscf.get('block') or []:
            block_data = self._parser_results(block)
            multiplicity = self._scalar(block_data.get('multiplicity'))
            weights = block_data.get('root_weights') or []
            if isinstance(weights, np.ndarray):
                weights = weights.tolist()
            n_roots = self._scalar(block_data.get('n_roots'))
            if multiplicity is not None:
                state_multiplicities.append(int(multiplicity))
                n_roots_per_multiplicity.append(len(weights) or int(n_roots or 0))
            state_weights.extend(float(weight) for weight in weights)
        return {
            'state_treatment': (
                'state_averaged' if len(state_weights) > 1 else 'state_specific'
            ),
            'n_state_groups': len(state_multiplicities) or None,
            'state_multiplicities': state_multiplicities or None,
            'n_roots_per_multiplicity': n_roots_per_multiplicity or None,
            'state_weights': state_weights or None,
        }

    def _is_casci(self, casscf: dict[str, Any], input_file: Any) -> bool:
        if isinstance(input_file, (list, tuple, np.ndarray)):
            input_file = ' '.join(str(item) for item in input_file)
        return 'MAXITER 1' in str(input_file or '').upper() or any(
            self._parser_results(block).get('casci_marker') is not None
            for block in casscf.get('block') or []
        )

    def _get_multireference_method_data(
        self,
        casscf: dict[str, Any],
        input_file: Any,
        method_type: str | None = None,
    ) -> dict[str, Any] | None:
        if not casscf:
            return None

        method = {
            'type': 'CASCI' if self._is_casci(casscf, input_file) else 'CASSCF',
            'active_space': self._build_active_space_data(casscf),
            **self._collect_state_data(casscf),
        }

        if method_type is not None and method['type'] != method_type:
            return None

        return method

    def get_multireference_scf_methods(
        self, casscf: dict[str, Any], input_file: Any
    ) -> list[dict[str, Any]]:
        method = self._get_multireference_method_data(casscf, input_file, 'CASSCF')
        if method is None:
            return []
        self._method = 'MultireferenceSCF'
        return [method]

    def get_multireference_ci_methods(
        self, casscf: dict[str, Any], input_file: Any
    ) -> list[dict[str, Any]]:
        method = self._get_multireference_method_data(casscf, input_file, 'CASCI')
        if method is None:
            return []
        self._method = 'MultireferenceCI'
        return [method]

    def get_multireference_pt_methods(
        self, casscf: dict[str, Any], input_file: Any
    ) -> list[dict[str, Any]]:
        method = self._get_multireference_method_data(casscf, input_file)
        if method is None:
            return []

        hints = '\n'.join(
            str(value)
            for value in (
                self._scalar(casscf.get('pt_method')),
                self._scalar(casscf.get('qd_nevpt_type')),
                input_file,
            )
            if value
        ).upper()
        if 'NEVPT2' not in hints:
            return []

        name = '-'.join(
            part
            for part, present in (
                ('QD', 'QD' in hints),
                ('SC', 'SC_NEVPT2' in hints or 'SC-NEVPT2' in hints),
                ('NEVPT2', True),
            )
            if present
        )
        self._method = 'MultireferencePT'
        return [{**method, 'type': 'NEVPT', 'order': 2, 'name': name}]

    @staticmethod
    def _build_local_correlation_dict(local_type: str) -> dict[str, Any]:
        return {
            'type': local_type,
            'spaces': [
                {
                    'space_kind': 'local_virtual_space',
                    'virtual_space_type': 'LNO' if local_type == 'LNO' else 'PNO',
                    'excitation_order': 2,
                }
            ],
        }

    @staticmethod
    def _excitation_orders(cc_type: str | None) -> list[int] | None:
        if not cc_type or not cc_type.startswith('CC'):
            return None
        return [
            order
            for marker, order in (('S', 1), ('D', 2), ('T', 3), ('Q', 4))
            if marker in cc_type[2:]
        ] or None

    @staticmethod
    def _infer_local_correlation_type(
        input_file: Any, cc_data: dict[str, Any]
    ) -> str | None:
        kc_formation = OutParser._scalar(cc_data.get('kc_formation')) or ''
        hints = f'{input_file or ""}\n{kc_formation}'.upper()
        return next(
            (
                local_type
                for local_type in ('DLPNO', 'LPNO', 'LNO', 'PNO')
                if local_type in hints
            ),
            None,
        )

    @staticmethod
    def _local_thresholds(
        block: dict[str, Any], include_triples: bool = True
    ) -> list[dict[str, Any]]:
        mapping = [
            ('tCutPairs', 'TCutPairs', 'pair_screening'),
            ('tCutPNO', 'TCutPNO', 'local_virtual_space'),
            ('tCutPNOSingles', 'TCutPNOSingles', 'local_virtual_space'),
            ('tCutMP2Pairs', 'TCutMP2Pairs', 'pair_screening'),
            ('tCutMKN', 'TCutMKN', 'occupied_domain'),
            ('tCutPAO', 'TCutPAO', 'local_virtual_space'),
            ('tCutEN', 'TCutEN', 'occupied_domain'),
            ('paoOverlapThresh', 'PAOOverlapThresh', 'local_virtual_space'),
        ]
        if include_triples:
            mapping.extend(
                [
                    ('tCutTNO', 'TCutTNO', 'local_virtual_space'),
                    ('tCutDOStrong', 'TCutDOStrong', 'occupied_domain'),
                    ('tCutMKNStrong', 'TCutMKNStrong', 'occupied_domain'),
                    ('tCutMKNWeak', 'TCutMKNWeak', 'occupied_domain'),
                    ('tCutDOWeak', 'TCutDOWeak', 'occupied_domain'),
                ]
            )
        return [
            {'name': name, 'value': OutParser._scalar(block[key]), 'applies_to': target}
            for key, name, target in mapping
            if block.get(key) is not None
        ]

    def get_orbital_localization_methods(
        self, loc_data: Any
    ) -> list[dict[str, Any]]:
        blocks = (
            loc_data if isinstance(loc_data, (list, tuple, np.ndarray)) else [loc_data]
        )
        methods = []
        for block in blocks:
            loc = self._parser_results(block)
            method = self._normalize_localization_method(loc.get('type'))
            if not method:
                continue
            method_data = {'method': method}
            orbital_range = loc.get('orbital_range')
            if isinstance(orbital_range, np.ndarray):
                orbital_range = orbital_range.tolist()
            if (
                isinstance(orbital_range, (list, tuple))
                and len(orbital_range) >= ORBITAL_RANGE_BOUNDS
            ):
                start, end = map(int, orbital_range[:ORBITAL_RANGE_BOUNDS])
                if end >= start:
                    method_data['n_localized_orbitals'] = end - start + 1
            methods.append(method_data)
        if methods:
            self._method = 'OrbitalLocalization'
        return sorted(
            methods, key=lambda item: item.get('n_localized_orbitals', 0), reverse=True
        )[:1]

    def get_perturbation_methods(
        self, mdci_data: dict[str, Any], mp2_data: dict[str, Any], input_file: Any
    ) -> list[dict[str, Any]]:
        mdci_data = mdci_data or {}
        if not mdci_data and not mp2_data:
            return []

        local_type = self._infer_local_correlation_type(input_file, mdci_data)
        method: dict[str, Any] = {'type': 'MP', 'order': 2}
        scaling = self._scalar(mdci_data.get('spin_component_scaling'))
        if isinstance(scaling, str) and scaling.strip().upper() in {'SCS', 'SOS'}:
            method['spin_component_scaling'] = scaling.strip().upper()
        if local_type:
            method['local_correlation'] = self._build_local_correlation_dict(local_type)
        thresholds = self._local_thresholds(mdci_data, include_triples=False)
        if thresholds:
            method['numerical_settings'] = [{'screening_thresholds': thresholds}]
        self._method = 'PerturbationMethod'
        return [method]

    def get_coupled_cluster_methods(
        self, cc_data: dict[str, Any], input_file: Any
    ) -> list[dict[str, Any]]:
        cc_data = cc_data or {}
        cc_type = self._scalar(cc_data.get('coupled_cluster_type'))
        if not isinstance(cc_type, str) or not cc_type:
            return []

        method: dict[str, Any] = {
            'type': cc_type,
            'excitation_order': self._excitation_orders(cc_type),
        }
        if (
            str(
                self._scalar(cc_data.get('perturbative_triple_excitations_on_off'))
            ).upper()
            == 'ON'
        ):
            method.update(
                perturbative_correction='(T)', perturbative_correction_order=[3]
            )
        if str(self._scalar(cc_data.get('f12_correction_on_off'))).upper() == 'ON':
            method['explicit_correlation'] = 'F12'
        local_type = self._infer_local_correlation_type(input_file, cc_data)
        if local_type:
            method['local_correlation'] = self._build_local_correlation_dict(local_type)
        thresholds = self._local_thresholds(cc_data)
        if thresholds:
            method['numerical_settings'] = [{'screening_thresholds': thresholds}]
        self._method = 'CC'
        return [{key: value for key, value in method.items() if value is not None}]

    def get_hf_methods(self, cc_data: dict[str, Any]) -> list[dict[str, Any]]:
        cc_data = cc_data or {}
        reference_form = self._hf_reference_form(
            cc_data.get('cc_reference_wavefunction')
        )
        if reference_form:
            self._method = 'HF'
        return [{'reference_form': reference_form}] if reference_form else []

    def get_relativity_model(self, source: dict[str, Any]) -> dict[str, Any]:
        relativistic = self._parser_results(source.get('relativistic_hamiltonian'))
        scf_settings = self._get_scf_settings(source)

        raw_method = self._scalar(
            relativistic.get('method') or scf_settings.get('scalar_relativistic_method')
        )
        if not isinstance(raw_method, str):
            return {}

        normalized = raw_method.upper()
        approximation = next(
            (
                approximation
                for marker, approximation in (
                    ('DOUGLAS-KROLL-HESS', 'DKH'),
                    ('DKH', 'DKH'),
                    ('ZORA', 'ZORA'),
                    ('FORA', 'FORA'),
                    ('IORA', 'IORA'),
                    ('X2C', 'X2C'),
                    ('BSS', 'BSS'),
                    ('NESC', 'NESC'),
                    ('PAULI', 'Pauli'),
                    ('SOMF', 'SOMF'),
                )
                if marker in normalized
            ),
            None,
        )
        if approximation is None:
            return {}

        model = {'level': 'scalar', 'approximation': approximation}
        dkh_order = self._scalar(
            relativistic.get('dkh_order') or scf_settings.get('dkh_order')
        )
        if approximation == 'DKH' and dkh_order is not None:
            model['dkh_order'] = int(dkh_order)
        # TODO implement support in mapping parser
        return model

    def get_relativity_models(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        method = self._method
        electronic_methods = (
            'DFT',
            'HF',
            'PerturbationMethod',
            'CC',
            'MultireferenceSCF',
            'MultireferenceCI',
            'MultireferencePT',
        )
        if method not in electronic_methods:
            return []
        model = self.get_relativity_model(source)
        model.setdefault('basis_set', {})
        return [model]

    def get_basis_set_components(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        source = source.get('basis_set', self.data)
        if not source:
            return []
        names = self._parser_results(source.get('basis_set_name'))
        totals = self._parser_results(source.get('basis_set_total'))
        allowed_keys = self._BASIS_SET_KEYS.get(self._method, frozenset())
        if not allowed_keys:
            return []
        components = []
        for key, role in self._BASIS_SET_ROLES.items():
            if key not in allowed_keys:
                continue
            name = self._scalar(names.get(key))
            if not isinstance(name, str) or not name.strip():
                continue
            component = {
                'source_key': key,
                'basis_set': name.strip(),
                'type': 'GTO',
                'role': role,
            }
            total = self._scalar(totals.get(key))
            if total is not None:
                component['n_total_basis_functions'] = int(total)
            components.append(component)

        return components

    def get_molecular_orbitals(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        self_consistent = self._navigate(src, 'single_point', 'self_consistent')
        basis_set_total = self._parser_results(src.get('basis_set_total'))

        orbital_energies = self._scalar(self_consistent.get('orbital_energies'))
        coefficients = str_to_mo_coefficients(
            self_consistent.get('molecular_orbital_coefficients')
        )
        if orbital_energies is None and coefficients is None:
            return []

        table = None
        if orbital_energies is not None:
            table = np.asarray(orbital_energies, dtype=np.float64)
            if (
                table.ndim != ORBITAL_ENERGY_N_DIMENSIONS
                or table.shape[1] < ORBITAL_ENERGY_N_COLUMNS
            ):
                table = None

        n_ao = self._scalar(basis_set_total.get('main_basis_set'))
        molecular_orbitals = {
            'n_mo': int(table.shape[0])
            if table is not None
            else int(coefficients.shape[0]),
            'n_ao': int(n_ao)
            if n_ao is not None
            else int(coefficients.shape[1])
            if coefficients is not None
            else int(table.shape[0]),
            'kind': 'canonical',
        }
        if table is not None:
            molecular_orbitals.update(
                {
                    'occupations': table[:, ORBITAL_OCCUPATION_COLUMN],
                    'energies': table[:, ORBITAL_ENERGY_EV_COLUMN] * ureg.electron_volt,
                }
            )
        if coefficients is not None:
            molecular_orbitals['coefficients'] = coefficients

        return [molecular_orbitals]

    def get_outputs(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        molecular_orbitals = self.get_molecular_orbitals(src)
        if not molecular_orbitals:
            return []

        return [
            {
                'model_system_ref': '/data/model_system/0',
                'molecular_orbitals': molecular_orbitals,
            }
        ]

    def build_workflow(
        self, archive: 'EntryArchive', logger: 'BoundLogger'
    ) -> SerialWorkflow | None:
        simulation = archive.data
        methods = simulation.model_method or []
        if not any(isinstance(method, HF) for method in methods) or not any(
            isinstance(method, CC) for method in methods
        ):
            return None

        tasks = [Task(name='HF')]
        if any(isinstance(method, OrbitalLocalization) for method in methods):
            tasks.append(Task(name='Orbital localization'))
        if any(isinstance(method, PerturbationMethod) for method in methods):
            tasks.append(Task(name='Local MP2'))
        outputs = (
            [Link(name='Outputs', section=simulation.outputs[-1])]
            if simulation.outputs
            else []
        )
        cc = next(method for method in methods if isinstance(method, CC))
        tasks.append(
            Task(name='Local CC' if cc.local_correlation else 'CC', outputs=outputs)
        )
        archive.workflow2 = SerialWorkflow(tasks=tasks)
        archive.workflow2.normalize(archive, logger)
        return archive.workflow2


class OrcaParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        reload(orca)

        reader = OutParser()
        reader.filepath = mainfile
        archive.data = Simulation(program=Program(name='ORCA'))
        metainfo_parser = MetainfoParser(data_object=archive.data)
        metainfo_parser.annotation_key = orca.OUT_KEY
        metainfo_parser.max_nested_level = 3

        try:
            reader.convert(metainfo_parser)
            reader.build_workflow(archive, logger)
        finally:
            metainfo_parser.close()
            reader.close()
