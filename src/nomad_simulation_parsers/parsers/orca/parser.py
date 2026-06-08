import re
from importlib import reload
from typing import TYPE_CHECKING, Any

import numpy as np
from nomad.datamodel.metainfo.workflow import Link, Task
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.units import ureg
from nomad_simulations.schema_packages.basis_set import (
    AtomCenteredBasisSet,
    BasisSetContainer,
)
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_method import (
    CC,
    DFT,
    HF,
    ActiveSpace,
    LocalCorrelation,
    LocalCorrelationSpace,
    MultireferenceCI,
    MultireferencePT,
    MultireferenceSCF,
    OrbitalLocalization,
    PerturbationMethod,
)
from nomad_simulations.schema_packages.numerical_settings import (
    LocalCorrelationSettings,
    LocalCorrelationThreshold,
)
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
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
        if not tokens[index].isdigit():
            index += 1
            continue

        mo_indices = []
        while index < len(tokens) and tokens[index].isdigit():
            mo_indices.append(int(tokens[index]))
            index += 1

        n_indices = len(mo_indices)
        index += n_indices * 2
        separator_tokens = tokens[index : index + n_indices]
        index += sum(1 for token in separator_tokens if token.startswith('-'))

        index = _read_coefficient_rows(tokens, index, mo_indices, coefficients_by_mo)

    return _coefficient_matrix(coefficients_by_mo)


def _read_coefficient_rows(
    tokens: list[str],
    index: int,
    mo_indices: list[int],
    coefficients_by_mo: dict[int, list[float]],
) -> int:
    n_indices = len(mo_indices)
    while index < len(tokens) and AO_ROW_RE.fullmatch(tokens[index]):
        row_start = index + 2
        row_end = row_start + n_indices
        if row_end > len(tokens):
            break

        try:
            row_values = [float(token) for token in tokens[row_start:row_end]]
        except ValueError:
            break

        for mo_index, coefficient in zip(mo_indices, row_values):
            coefficients_by_mo.setdefault(mo_index, []).append(coefficient)
        index = row_end

    return index


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
    def __init__(self) -> None:
        super().__init__(text_parser=OutReader())
        self.parse_only_required = False
        self.text_parser.parse_only_required = False
        self.text_parser.findlazy = False

    def load_file(self) -> OutReader:
        if self.filepath:
            self.text_parser.findlazy = False
            self.text_parser.mainfile = self.filepath
        return self.text_parser

    @staticmethod
    def _as_dict(value: Any) -> dict[str, Any]:
        return value._results if hasattr(value, '_results') else value or {}

    def _get_mdci_data(self, source: dict[str, Any]) -> dict[str, Any]:
        single_point = self._as_dict(source.get('single_point'))
        return self._as_dict(single_point.get('ci'))

    @staticmethod
    def _scalar(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(value, (list, tuple)):
            return value[0] if value else None
        return value

    def get_program_name(self, src: dict[str, Any]) -> str:
        return 'ORCA'

    def get_atoms(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        single_point = self._as_dict(src.get('single_point', {}))
        coordinates = single_point.get('cartesian_coordinates', [])
        if not coordinates:
            return []

        symbols, positions = str_to_cartesian_coordinates(coordinates)
        if not symbols:
            return []

        system = {
            'is_representative': True,
            'positions': positions,
            'particle_states': [{'chemical_symbol': symbol} for symbol in symbols],
        }

        self_consistent = self._as_dict(single_point.get('self_consistent', {}))
        scf_settings = self._as_dict(self_consistent.get('scf_settings', {}))
        total_charge = self._scalar(scf_settings.get('total_charge'))
        if total_charge is not None:
            system['total_charge'] = int(total_charge)

        multiplicity = self._scalar(scf_settings.get('multiplicity'))
        if multiplicity is not None:
            system['total_spin'] = int(multiplicity) - 1

        return [system]

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

    def get_dft(self, src: dict[str, Any]) -> dict[str, Any]:
        single_point = self._as_dict(src.get('single_point'))
        self_consistent = self._as_dict(single_point.get('self_consistent'))
        scf_settings = self._as_dict(self_consistent.get('scf_settings'))
        if not scf_settings:
            return {}

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

        xc = {}
        if functional_key:
            xc['functional_key'] = functional_key
        hf_fraction = self._scalar(scf_settings.get('fraction_hf_exchange'))
        if hf_fraction is not None:
            xc['global_exact_exchange'] = hf_fraction
        if not xc:
            return {}

        method = {'xc': xc}
        reference_form = self._dft_reference_form(scf_settings.get('hf_type'))
        if reference_form:
            method['reference_form'] = reference_form
        return method

    def _get_multireference_method_data(self, source: dict[str, Any]) -> dict[str, Any]:
        single_point = self._as_dict(source.get('single_point'))
        casscf = self._as_dict(single_point.get('casscf'))
        if not casscf:
            return {}

        active_space = {
            'n_active_electrons': self._scalar(casscf.get('n_active_electrons')),
            'n_active_orbitals': self._scalar(casscf.get('n_active_orbitals')),
            'orbital_space_type': 'CAS',
        }
        state_multiplicities = []
        n_roots_per_multiplicity = []
        state_weights = []
        for block in casscf.get('block') or []:
            block_data = self._as_dict(block)
            multiplicity = self._scalar(block_data.get('multiplicity'))
            weights = block_data.get('root_weights') or []
            if isinstance(weights, np.ndarray):
                weights = weights.tolist()
            n_roots = self._scalar(block_data.get('n_roots'))
            if multiplicity is not None:
                state_multiplicities.append(int(multiplicity))
                n_roots_per_multiplicity.append(len(weights) or int(n_roots or 0))
            state_weights.extend(float(weight) for weight in weights)

        input_file = source.get('input_file') or ''
        if isinstance(input_file, (list, tuple, np.ndarray)):
            input_file = ' '.join(str(item) for item in input_file)
        is_casci = 'MAXITER 1' in str(input_file).upper() or any(
            self._as_dict(block).get('casci_marker') is not None
            for block in casscf.get('block') or []
        )
        return {
            'type': 'CASCI' if is_casci else 'CASSCF',
            'active_space': active_space,
            'state_treatment': (
                'state_averaged' if len(state_weights) > 1 else 'state_specific'
            ),
            'n_state_groups': len(state_multiplicities) or None,
            'state_multiplicities': state_multiplicities or None,
            'n_roots_per_multiplicity': n_roots_per_multiplicity or None,
            'state_weights': state_weights or None,
        }

    def get_multireference_scf_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        method = self._get_multireference_method_data(source)
        return [method] if method.get('type') == 'CASSCF' else []

    def get_multireference_ci_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        method = self._get_multireference_method_data(source)
        return [method] if method.get('type') == 'CASCI' else []

    def get_multireference_pt_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        method = self._get_multireference_method_data(source)
        if not method:
            return []

        casscf = self._as_dict(self._as_dict(source.get('single_point')).get('casscf'))
        input_file = source.get('input_file') or ''
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
        return [{**method, 'type': 'NEVPT', 'order': 2, 'name': name}]

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
        source: dict[str, Any], cc_data: dict[str, Any]
    ) -> str | None:
        input_file = source.get('input_file') or ''
        kc_formation = OutParser._scalar(cc_data.get('kc_formation')) or ''
        hints = f'{input_file}\n{kc_formation}'.upper()
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
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        loc_data = self._as_dict(source.get('single_point')).get('loc')
        blocks = (
            loc_data if isinstance(loc_data, (list, tuple, np.ndarray)) else [loc_data]
        )
        methods = []
        for block in blocks:
            loc = self._as_dict(block)
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
        return sorted(
            methods, key=lambda item: item.get('n_localized_orbitals', 0), reverse=True
        )[:1]

    def get_perturbation_methods(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        single_point = self._as_dict(source.get('single_point'))
        mdci_data = self._get_mdci_data(source)
        mp2_data = self._as_dict(single_point.get('mp2'))
        if not mdci_data and not mp2_data:
            return []

        local_type = self._infer_local_correlation_type(source, mdci_data)
        method: dict[str, Any] = {'type': 'MP', 'order': 2}
        scaling = self._scalar(mdci_data.get('spin_component_scaling'))
        if isinstance(scaling, str) and scaling.strip().upper() in {'SCS', 'SOS'}:
            method['spin_component_scaling'] = scaling.strip().upper()
        if local_type:
            method['local_correlation'] = {
                'type': local_type,
                'spaces': [
                    {
                        'space_kind': 'local_virtual_space',
                        'virtual_space_type': 'LNO' if local_type == 'LNO' else 'PNO',
                        'excitation_order': 2,
                    }
                ],
            }
        thresholds = self._local_thresholds(mdci_data, include_triples=False)
        if thresholds:
            method['numerical_settings'] = [{'screening_thresholds': thresholds}]
        return [method]

    def get_coupled_cluster_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        cc_data = self._get_mdci_data(source)
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
        local_type = self._infer_local_correlation_type(source, cc_data)
        if local_type:
            method['local_correlation'] = {
                'type': local_type,
                'spaces': [
                    {
                        'space_kind': 'local_virtual_space',
                        'virtual_space_type': 'LNO' if local_type == 'LNO' else 'PNO',
                        'excitation_order': 2,
                    }
                ],
            }
        thresholds = self._local_thresholds(cc_data)
        if thresholds:
            method['numerical_settings'] = [{'screening_thresholds': thresholds}]
        return [{key: value for key, value in method.items() if value is not None}]

    def get_hf_methods(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        cc_data = self._get_mdci_data(source)
        reference_form = self._hf_reference_form(
            cc_data.get('cc_reference_wavefunction')
        )
        return [{'reference_form': reference_form}] if reference_form else []

    def get_basis_set_components(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        names = self._as_dict(source.get('basis_set_name'))
        totals = self._as_dict(source.get('basis_set_total'))
        roles = {
            'main_basis_set': 'orbital',
            'auxj_basis_set': 'auxiliary_scf',
            'auxjk_basis_set': 'auxiliary_scf',
            'auxc_basis_set': 'auxiliary_post_hf',
        }
        components = []
        for key, role in roles.items():
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
        single_point = self._as_dict(src.get('single_point', {}))
        self_consistent = self._as_dict(single_point.get('self_consistent', {}))
        basis_set_total = self._as_dict(src.get('basis_set_total', {}))

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
            'm_def': MolecularOrbitals.m_def.qualified_name(),
            'n_mo': int(table.shape[0])
            if table is not None
            else int(coefficients.shape[0]),
            'n_ao': int(n_ao)
            if n_ao is not None
            else int(coefficients.shape[1])
            if coefficients is not None
            else int(table.shape[0]),
            'mo_type': 'canonical',
        }
        if table is not None:
            molecular_orbitals.update(
                {
                    'mo_occupations': table[:, ORBITAL_OCCUPATION_COLUMN],
                    'mo_energies': table[:, ORBITAL_ENERGY_EV_COLUMN]
                    * ureg.electron_volt,
                }
            )
        if coefficients is not None:
            molecular_orbitals['mo_coefficients'] = coefficients

        return [molecular_orbitals]

    def get_outputs(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        molecular_orbitals = self.get_molecular_orbitals(src)
        if not molecular_orbitals:
            return []

        return [
            {
                'model_system_ref': '/data/model_system/0',
                'electronic_eigenvalues': molecular_orbitals,
            }
        ]

    @staticmethod
    def _build_active_space(data: dict[str, Any] | None) -> ActiveSpace | None:
        return ActiveSpace(**data) if data else None

    @staticmethod
    def _build_local_correlation(data: dict[str, Any]) -> LocalCorrelation:
        return LocalCorrelation(
            type=data.get('type'),
            spaces=[
                LocalCorrelationSpace(**space_data)
                for space_data in data.get('spaces', [])
            ],
        )

    @staticmethod
    def _build_local_settings(data: dict[str, Any]) -> LocalCorrelationSettings:
        return LocalCorrelationSettings(
            screening_thresholds=[
                LocalCorrelationThreshold(**threshold)
                for threshold in data.get('screening_thresholds', [])
            ]
        )

    @classmethod
    def _build_perturbation(cls, data: dict[str, Any]) -> PerturbationMethod:
        return PerturbationMethod(
            type=data.get('type'),
            order=data.get('order'),
            spin_component_scaling=data.get('spin_component_scaling'),
            local_correlation=cls._build_local_correlation(data['local_correlation'])
            if data.get('local_correlation')
            else None,
            numerical_settings=[
                cls._build_local_settings(settings)
                for settings in data.get('numerical_settings', [])
            ],
        )

    @classmethod
    def _build_cc(cls, data: dict[str, Any]) -> CC:
        return CC(
            type=data.get('type'),
            excitation_order=data.get('excitation_order'),
            perturbative_correction=data.get('perturbative_correction'),
            perturbative_correction_order=data.get('perturbative_correction_order'),
            explicit_correlation=data.get('explicit_correlation'),
            local_correlation=cls._build_local_correlation(data['local_correlation'])
            if data.get('local_correlation')
            else None,
            numerical_settings=[
                cls._build_local_settings(settings)
                for settings in data.get('numerical_settings', [])
            ],
        )

    def enrich_methods(self, simulation: Simulation) -> None:
        source = self.text_parser.results or {}
        existing = list(simulation.model_method or [])
        additions = []

        dft_data = self.get_dft(source)
        if dft_data:
            additions.append(DFT(**dft_data))

        for data in self.get_hf_methods(source):
            additions.append(HF(**data))

        localization_data = self.get_orbital_localization_methods(source)
        localization = (
            OrbitalLocalization(**localization_data[0]) if localization_data else None
        )
        if localization:
            additions.append(localization)

        for data in self.get_perturbation_methods(source):
            additions.append(self._build_perturbation(data))

        for data in self.get_coupled_cluster_methods(source):
            additions.append(self._build_cc(data))

        multireference = []
        for data in self.get_multireference_scf_methods(source):
            multireference.append(
                MultireferenceSCF(
                    **{
                        **data,
                        'active_space': self._build_active_space(
                            data.get('active_space')
                        ),
                    }
                )
            )
        for data in self.get_multireference_ci_methods(source):
            multireference.append(
                MultireferenceCI(
                    **{
                        **data,
                        'active_space': self._build_active_space(
                            data.get('active_space')
                        ),
                    }
                )
            )
        for data in self.get_multireference_pt_methods(source):
            multireference.append(
                MultireferencePT(
                    **{
                        **data,
                        'active_space': self._build_active_space(
                            data.get('active_space')
                        ),
                    }
                )
            )

        replaced_types = (
            DFT,
            HF,
            OrbitalLocalization,
            PerturbationMethod,
            CC,
            MultireferenceSCF,
            MultireferenceCI,
            MultireferencePT,
        )
        preserved = [
            method for method in existing if not isinstance(method, replaced_types)
        ]
        simulation.model_method = additions + preserved + multireference

    def enrich_basis_sets(self, simulation: Simulation) -> None:
        components = self.get_basis_set_components(self.text_parser.results or {})
        if not components:
            return

        species_scope = (
            list(simulation.model_system[0].particle_states)
            if simulation.model_system and simulation.model_system[0].particle_states
            else None
        )
        component_lookup = {
            component['source_key']: component for component in components
        }
        for method in simulation.model_method or []:
            if isinstance(method, (HF, DFT)):
                keys = ('main_basis_set', 'auxj_basis_set', 'auxjk_basis_set')
            elif isinstance(method, (PerturbationMethod, CC)):
                keys = ('main_basis_set', 'auxc_basis_set')
            else:
                continue

            selected = [
                component_lookup[key] for key in keys if key in component_lookup
            ]
            if not selected:
                continue
            basis_components = []
            for component in selected:
                data = {
                    key: value
                    for key, value in component.items()
                    if key != 'source_key'
                }
                if species_scope:
                    data['species_scope'] = species_scope
                basis_components.append(AtomCenteredBasisSet(**data))
            method.numerical_settings.append(
                BasisSetContainer(basis_set_components=basis_components)
            )

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
        archive.data = Simulation()
        metainfo_parser = MetainfoParser(data_object=archive.data)
        metainfo_parser.annotation_key = orca.OUT_KEY
        metainfo_parser.max_nested_level = 3

        try:
            reader.convert(metainfo_parser)
            reader.enrich_methods(archive.data)
            reader.enrich_basis_sets(archive.data)
            reader.build_workflow(archive, logger)
        finally:
            metainfo_parser.close()
            reader.close()
