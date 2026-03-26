from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from importlib import reload

import numpy as np
from nomad.datamodel.metainfo.workflow import Link, Task
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.basis_set import (
    AtomCenteredBasisSet,
    BasisSetContainer,
)
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_method import (
    CC,
    DFT,
    HF,
    LocalCorrelation,
    LocalCorrelationSpace,
    OrbitalLocalization,
    PerturbationMethod,
)
from nomad_simulations.schema_packages.numerical_settings import (
    LocalCorrelationSettings,
    LocalCorrelationThreshold,
)
from nomad_simulations.schema_packages.workflow.general import SerialWorkflow

from nomad_simulation_parsers.schema_packages import orca
from nomad_simulation_parsers.schema_packages.utils import remove_mapping_annotations

from .text_parser import OutReader

LOGGER = get_logger(__name__)


def str_to_cartesian_coordinates(val_in):
    val_in_cleaned = [
        val.replace('>', '') if isinstance(val, str) else val
        for val in val_in
        if val != '>'
    ]

    if isinstance(val_in_cleaned, list):
        symbols = []
        coordinates = []
        for i in range(0, len(val_in_cleaned), 4):
            symbol = val_in_cleaned[i]
            if isinstance(symbol, str):
                symbol = symbol.replace('>', '')
            symbols.append(symbol)
            coordinates.append(val_in_cleaned[i + 1 : i + 4])
            # print(coordinates)
        coordinates = np.array(coordinates, dtype=float)
        return symbols, coordinates


class OutParser(MappingTextParser):
    """
    Couples OrcaTextParser (regex) with a few convenience getters that
    the mapping rules will call.
    """

    def __init__(self):
        super().__init__(text_parser=OutReader())
        # Parse all quantities so helper-driven mappings (e.g., CASSCF) are available
        self.parse_only_required = False
        self.text_parser.parse_only_required = False
        self.text_parser.findlazy = False

    def load_file(self):
        if self.filepath:
            self.text_parser.findlazy = False
            self.text_parser.mainfile = self.filepath
        return self.text_parser

    def get_program_data(self, src: dict[str, Any]) -> dict[str, Any]:
        return {
            'program_name': 'ORCA',
            'program_version': src.get('program_version'),
        }

    def get_atoms(self, src: dict[str, Any]):
        # ← revert to the original nested lookup
        coords = src.get('single_point', {}).get('cartesian_coordinates', [])
        if not coords:
            return []

        syms, pos = str_to_cartesian_coordinates(coords)
        atoms = [{'chemical_symbol': s} for s in syms]
        return [{'positions': pos, 'particle_states': atoms}]

    def get_dft(self, src: dict[str, Any]) -> dict[str, Any]:
        """
        Build DFT model method data from ORCA SCF settings.
        """
        scf_settings = (
            src.get('single_point', {})
            .get('self_consistent', {})
            .get('scf_settings', {})
        )
        if not scf_settings:
            return {}

        exchange = scf_settings.get('exchange_functional')
        correlation = scf_settings.get('correlation_functional') or scf_settings.get(
            'correl_functional'
        )
        hf_frac = scf_settings.get('fraction_hf_exchange')

        functional_key = None
        if isinstance(exchange, str) and isinstance(correlation, str):
            if exchange and correlation:
                functional_key = (
                    exchange if exchange == correlation else f'{exchange}+{correlation}'
                )
        elif isinstance(exchange, str) and exchange:
            functional_key = exchange
        elif isinstance(correlation, str) and correlation:
            functional_key = correlation

        xc = {}
        if functional_key:
            xc['functional_key'] = functional_key
        if hf_frac is not None:
            xc['global_exact_exchange'] = hf_frac
        if not xc:
            return {}

        return {'xc': xc}

    def get_numerical_settings(self, source: dict[str, Any]) -> dict[str, Any]:
        scf_convergence = (
            source.get('single_point', {})
            .get('self_consistent', {})
            .get('scf_settings', {})
        )

        return {
            'n_max_iterations': scf_convergence.get('n_max_iterations', 2575),
            'threshold_change': scf_convergence.get('energy_change_tolerance', 1e-8),
        }

    def get_multireference_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        casscf = source.get('single_point', {}).get('casscf') if source else None
        if hasattr(casscf, '_results'):
            casscf = casscf._results
        if not casscf:
            return []

        active_space = {
            'n_active_electrons': casscf.get('n_active_electrons'),
            'n_active_orbitals': casscf.get('n_active_orbitals'),
            'orbital_space_type': 'CAS',
        }
        active_space = {k: v for k, v in active_space.items() if v is not None}

        state_multiplicities: list[int] = []
        n_roots_per_multiplicity: list[int] = []
        state_weights: list[float] = []
        for block in casscf.get('block') or []:
            block_data = block._results if hasattr(block, '_results') else block
            multiplicity = block_data.get('multiplicity')
            weights = block_data.get('root_weights') or []
            n_roots = block_data.get('n_roots')

            state_multiplicities.append(multiplicity)
            n_roots_per_multiplicity.append(len(weights) or n_roots)
            state_weights.extend(weights)

        reference_type = (
            'state_averaged'
            if state_weights and len(state_weights) > 1
            else 'state_specific'
        )
        n_state_groups = len(state_multiplicities) if state_multiplicities else None

        return [
            {
                'type': 'CASSCF',
                'active_space': active_space or None,
                'reference_type': reference_type,
                'n_state_groups': n_state_groups,
                'state_multiplicities': state_multiplicities or None,
                'n_roots_per_multiplicity': n_roots_per_multiplicity or None,
                'state_weights': state_weights or None,
            }
        ]

    def get_basis_set_components(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        basis_set_names = self._as_dict(source.get('basis_set_name'))
        basis_set_totals = self._as_dict(source.get('basis_set_total'))
        if not basis_set_names:
            return []

        role_mapping = {
            'main_basis_set': 'orbital',
            'auxj_basis_set': 'auxiliary_scf',
            'auxjk_basis_set': 'auxiliary_scf',
            'auxc_basis_set': 'auxiliary_post_hf',
        }

        basis_components = []
        for key, role in role_mapping.items():
            basis_name = self._scalar(basis_set_names.get(key))
            if not isinstance(basis_name, str) or not basis_name.strip():
                continue

            component = {
                'source_key': key,
                'basis_set': basis_name.strip(),
                'type': 'GTO',
                'role': role,
            }

            n_total_basis_functions = self._scalar(basis_set_totals.get(key))
            if n_total_basis_functions is not None:
                component['n_total_basis_functions'] = int(n_total_basis_functions)

            basis_components.append(component)

        return basis_components

    @staticmethod
    def _as_dict(value: Any) -> dict[str, Any]:
        return value._results if hasattr(value, '_results') else value or {}

    @staticmethod
    def _scalar(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(value, (list, tuple)):
            return value[0] if value else None
        return value

    @staticmethod
    def _normalize_localization_method(value: Any) -> str | None:
        raw = OutParser._scalar(value)
        if not isinstance(raw, str):
            return None

        raw_upper = raw.upper()
        mapping = {
            'FOSTER-BOYS': 'Foster-Boys',
            'PIPEK-MEZEY': 'Pipek-Mezey',
            'EDMISTON-RUEDENBERG': 'Edmiston-Ruedenberg',
            'IBO': 'IBO',
            'AIOPM-NBO': 'AIOPM-NBO',
            'NBO': 'NBO',
        }
        for key, normalized in mapping.items():
            if key in raw_upper:
                return normalized
        return None

    @staticmethod
    def _determinant_from_reference(reference_wavefunction: Any) -> str | None:
        reference = OutParser._scalar(reference_wavefunction)
        if not isinstance(reference, str):
            return None

        reference = reference.upper()
        if reference == 'RHF':
            return 'restricted'
        if reference == 'UHF':
            return 'unrestricted'
        if reference == 'ROHF':
            return 'restricted-open-shell'
        return None

    @staticmethod
    def _excitation_orders_from_cc_type(cc_type: str | None) -> list[int] | None:
        if not cc_type or not cc_type.startswith('CC'):
            return None

        suffix = cc_type[2:]
        orders: list[int] = []
        if 'S' in suffix:
            orders.append(1)
        if 'D' in suffix:
            orders.append(2)
        if 'T' in suffix:
            orders.append(3)
        if 'Q' in suffix:
            orders.append(4)
        return orders or None

    @staticmethod
    def _infer_local_correlation_type(
        source: dict[str, Any], cc_data: dict[str, Any]
    ) -> str | None:
        input_file = source.get('input_file') or ''
        kc_formation = OutParser._scalar(cc_data.get('kc_formation')) or ''
        hints = f'{input_file}\n{kc_formation}'.upper()

        for local_type in ('DLPNO', 'LPNO', 'LNO', 'PNO'):
            if local_type in hints:
                return local_type
        return None

    @staticmethod
    def _virtual_space_type_from_local_type(local_type: str | None) -> str | None:
        if local_type == 'LNO':
            return 'LNO'
        if local_type in {'DLPNO', 'LPNO', 'PNO'}:
            return 'PNO'
        return None

    @staticmethod
    def _normalize_spin_component_scaling(value: Any) -> str | None:
        raw = OutParser._scalar(value)
        if not isinstance(raw, str):
            return None

        normalized = raw.strip().upper()
        if normalized == 'SCS':
            return 'SCS'
        if normalized == 'SOS':
            return 'SOS'
        return None

    @staticmethod
    def _local_thresholds(
        block_data: dict[str, Any],
        threshold_mapping: tuple[tuple[str, str, str], ...],
    ) -> list[dict[str, Any]]:
        thresholds = []
        for source_key, threshold_name, applies_to in threshold_mapping:
            value = OutParser._scalar(block_data.get(source_key))
            if value is None:
                continue
            thresholds.append(
                {
                    'name': threshold_name,
                    'value': value,
                    'applies_to': applies_to,
                }
            )
        return thresholds

    def get_orbital_localization_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        loc_data = source.get('single_point', {}).get('loc')
        if not loc_data:
            return []

        loc_blocks = (
            list(loc_data)
            if isinstance(loc_data, (list, tuple, np.ndarray))
            else [loc_data]
        )
        methods = []
        for loc_block in loc_blocks:
            loc = self._as_dict(loc_block)
            if not loc:
                continue

            method = self._normalize_localization_method(loc.get('type'))
            if method is None:
                continue

            orbital_range = loc.get('orbital_range')
            n_localized_orbitals = None
            if isinstance(orbital_range, np.ndarray):
                orbital_range = orbital_range.tolist()
            if isinstance(orbital_range, (list, tuple)) and len(orbital_range) >= 2:
                try:
                    start = int(orbital_range[0])
                    end = int(orbital_range[1])
                except (TypeError, ValueError):
                    start = end = None
                if start is not None and end is not None and end >= start:
                    n_localized_orbitals = end - start + 1

            method_data = {'method': method}
            if n_localized_orbitals is not None:
                method_data['n_localized_orbitals'] = n_localized_orbitals
            methods.append(method_data)

        if not methods:
            return []

        methods.sort(
            key=lambda method_data: method_data.get('n_localized_orbitals') or 0,
            reverse=True,
        )
        return [methods[0]]

    def get_perturbation_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        cc_data = self._as_dict(source.get('single_point', {}).get('cc'))
        ci_data = self._as_dict(source.get('single_point', {}).get('ci'))
        if not ci_data:
            return []

        if ci_data.get('mp2_total_energy') is None and ci_data.get(
            'sl_mp2_correlation_energy'
        ) is None:
            return []

        local_type = self._infer_local_correlation_type(source, cc_data)
        method: dict[str, Any] = {
            'type': 'MP',
            'order': 2,
            'determinant': self._determinant_from_reference(
                cc_data.get('cc_reference_wavefunction')
            ),
            'spin_component_scaling': self._normalize_spin_component_scaling(
                ci_data.get('spin_component_scaling')
            ),
        }

        if local_type:
            virtual_space_type = self._virtual_space_type_from_local_type(local_type)
            method['local_correlation'] = {
                'type': local_type,
                'spaces': [
                    {
                        'kind': 'orbital',
                        'virtual_space_type': virtual_space_type,
                        'excitation_order': 2,
                    }
                ],
            }

        thresholds = []
        thresholds.extend(
            self._local_thresholds(
                cc_data,
                (
                    ('tCutPairs', 'TCutPairs', 'pair_screening'),
                    ('tCutPNO', 'TCutPNO', 'orbital'),
                    ('tCutPNOSingles', 'TCutPNOSingles', 'orbital'),
                ),
            )
        )
        thresholds.extend(
            self._local_thresholds(
                ci_data,
                (('tCutMP2Pairs', 'TCutMP2Pairs', 'pair_screening'),),
            )
        )
        if thresholds:
            method['numerical_settings'] = [{'screening_thresholds': thresholds}]

        method = {key: value for key, value in method.items() if value is not None}
        return [method]

    def get_coupled_cluster_methods(
        self, source: dict[str, Any]
    ) -> list[dict[str, Any]]:
        cc_data = self._as_dict(source.get('single_point', {}).get('cc'))
        if not cc_data:
            return []

        cc_type = self._scalar(cc_data.get('coupled_cluster_type'))
        if not isinstance(cc_type, str) or not cc_type:
            return []

        method: dict[str, Any] = {
            'type': cc_type,
            'excitation_order': self._excitation_orders_from_cc_type(cc_type),
            'determinant': self._determinant_from_reference(
                cc_data.get('cc_reference_wavefunction')
            ),
        }

        perturbative_triples = self._scalar(
            cc_data.get('perturbative_triple_excitations_on_off')
        )
        if isinstance(perturbative_triples, str) and perturbative_triples.upper() == 'ON':
            method['perturbative_correction'] = '(T)'
            method['perturbative_correction_order'] = [3]

        f12_correction = self._scalar(cc_data.get('f12_correction_on_off'))
        if isinstance(f12_correction, str) and f12_correction.upper() == 'ON':
            method['explicit_correlation'] = 'F12'

        local_type = self._infer_local_correlation_type(source, cc_data)
        if local_type:
            virtual_space_type = self._virtual_space_type_from_local_type(local_type)
            method['local_correlation'] = {
                'type': local_type,
                'spaces': [
                    {
                        'kind': 'orbital',
                        'virtual_space_type': virtual_space_type,
                        'excitation_order': 2,
                    }
                ],
            }

        thresholds = self._local_thresholds(
            cc_data,
            (
                ('tCutPairs', 'TCutPairs', 'pair_screening'),
                ('tCutPNO', 'TCutPNO', 'orbital'),
                ('tCutPNOSingles', 'TCutPNOSingles', 'orbital'),
            ),
        )
        if thresholds:
            method['numerical_settings'] = [{'screening_thresholds': thresholds}]

        method = {key: value for key, value in method.items() if value is not None}
        return [method]

    def get_hf_methods(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        cc_data = self._as_dict(source.get('single_point', {}).get('cc'))
        if not cc_data:
            return []

        reference = self._scalar(cc_data.get('cc_reference_wavefunction'))
        if reference not in {'RHF', 'UHF', 'ROHF'}:
            return []

        return [{'type': reference}]

    @staticmethod
    def _build_orbital_localization_section(
        method_data: dict[str, Any],
    ) -> OrbitalLocalization:
        return OrbitalLocalization(**method_data)

    @staticmethod
    def _build_hf_section(method_data: dict[str, Any]) -> HF:
        return HF(type=method_data.get('type'))

    @staticmethod
    def _build_dft_section(method_data: dict[str, Any]) -> DFT:
        return DFT(**method_data)

    @staticmethod
    def _build_perturbation_section(method_data: dict[str, Any]) -> PerturbationMethod:
        kwargs = {
            'type': method_data.get('type'),
            'order': method_data.get('order'),
            'determinant': method_data.get('determinant'),
            'spin_component_scaling': method_data.get('spin_component_scaling'),
        }

        local_correlation = method_data.get('local_correlation')
        if local_correlation:
            kwargs['local_correlation'] = OutParser._build_local_correlation_section(
                local_correlation
            )

        numerical_settings = method_data.get('numerical_settings') or []
        if numerical_settings:
            kwargs['numerical_settings'] = [
                OutParser._build_local_correlation_settings(settings_data)
                for settings_data in numerical_settings
            ]

        return PerturbationMethod(**kwargs)

    @staticmethod
    def _build_local_correlation_section(
        local_data: dict[str, Any],
    ) -> LocalCorrelation:
        spaces = [
            LocalCorrelationSpace(**space_data)
            for space_data in local_data.get('spaces', [])
        ]
        return LocalCorrelation(
            type=local_data.get('type'),
            spaces=spaces or None,
        )

    @staticmethod
    def _build_local_correlation_settings(
        settings_data: dict[str, Any],
    ) -> LocalCorrelationSettings:
        thresholds = [
            LocalCorrelationThreshold(**threshold_data)
            for threshold_data in settings_data.get('screening_thresholds', [])
        ]
        return LocalCorrelationSettings(screening_thresholds=thresholds or None)

    @staticmethod
    def _build_basis_set_container(
        component_data: list[dict[str, Any]], species_scope: list[Any] | None
    ) -> BasisSetContainer:
        basis_set_components = []
        for component in component_data:
            kwargs = {
                'basis_set': component.get('basis_set'),
                'type': component.get('type'),
                'role': component.get('role'),
                'n_total_basis_functions': component.get('n_total_basis_functions'),
            }
            if species_scope:
                kwargs['species_scope'] = species_scope
            basis_set_components.append(AtomCenteredBasisSet(**kwargs))

        return BasisSetContainer(basis_set_components=basis_set_components)

    def _build_cc_section(self, method_data: dict[str, Any]) -> CC:
        kwargs = {
            'type': method_data.get('type'),
            'determinant': method_data.get('determinant'),
            'excitation_order': method_data.get('excitation_order'),
            'perturbative_correction': method_data.get('perturbative_correction'),
            'perturbative_correction_order': method_data.get(
                'perturbative_correction_order'
            ),
            'explicit_correlation': method_data.get('explicit_correlation'),
        }

        local_correlation = method_data.get('local_correlation')
        if local_correlation:
            kwargs['local_correlation'] = self._build_local_correlation_section(
                local_correlation
            )

        numerical_settings = method_data.get('numerical_settings') or []
        if numerical_settings:
            kwargs['numerical_settings'] = [
                self._build_local_correlation_settings(settings_data)
                for settings_data in numerical_settings
            ]

        return CC(**kwargs)

    def enrich_local_correlation_methods(self, simulation: Simulation) -> None:
        if simulation is None:
            return

        source = self.text_parser.results or {}
        existing_methods = [
            method
            for method in (simulation.model_method or [])
            if not isinstance(method, (HF, OrbitalLocalization, PerturbationMethod, CC))
        ]
        dft_data = self.get_dft(source)
        if dft_data and not any(isinstance(method, DFT) for method in existing_methods):
            existing_methods.insert(0, self._build_dft_section(dft_data))

        hf_data = self.get_hf_methods(source)
        localization_data = self.get_orbital_localization_methods(source)
        mp2_data = self.get_perturbation_methods(source)
        cc_data = self.get_coupled_cluster_methods(source)
        method_sequence = []

        hf_method = self._build_hf_section(hf_data[0]) if hf_data else None
        if hf_method is not None:
            method_sequence.append(hf_method)

        localization = (
            self._build_orbital_localization_section(localization_data[0])
            if localization_data
            else None
        )
        if localization is not None:
            method_sequence.append(localization)

        mp2_method = (
            self._build_perturbation_section(mp2_data[0]) if mp2_data else None
        )
        if (
            mp2_method is not None
            and mp2_method.local_correlation is not None
            and localization is not None
        ):
            mp2_method.local_correlation.orbital_localization_ref = localization
        if mp2_method is not None:
            method_sequence.append(mp2_method)

        cc_method = self._build_cc_section(cc_data[0]) if cc_data else None
        if (
            cc_method is not None
            and cc_method.local_correlation is not None
            and localization is not None
        ):
            cc_method.local_correlation.orbital_localization_ref = localization
        if cc_method is not None:
            method_sequence.append(cc_method)

        method_sequence.extend(existing_methods)

        if method_sequence:
            simulation.model_method = method_sequence

    def enrich_basis_sets(self, simulation: Simulation) -> None:
        if simulation is None or not simulation.model_method:
            return

        basis_components = self.get_basis_set_components(self.text_parser.results or {})
        if not basis_components:
            return

        species_scope = None
        if simulation.model_system and simulation.model_system[0].particle_states:
            species_scope = list(simulation.model_system[0].particle_states)

        component_lookup = {
            component['source_key']: component for component in basis_components
        }

        for method in simulation.model_method:
            selected_keys = []
            if isinstance(method, (HF, DFT)):
                selected_keys = ['main_basis_set', 'auxj_basis_set', 'auxjk_basis_set']
            elif isinstance(method, (PerturbationMethod, CC)):
                selected_keys = ['main_basis_set', 'auxc_basis_set']

            selected_components = [
                component_lookup[key]
                for key in selected_keys
                if key in component_lookup
            ]
            if not selected_components:
                continue

            basis_container = self._build_basis_set_container(
                selected_components, species_scope
            )
            if method.numerical_settings is None:
                method.numerical_settings = []
            method.numerical_settings.append(basis_container)

    def build_workflow(
        self, simulation: Simulation, logger: 'BoundLogger'
    ):
        if simulation is None or not simulation.model_method:
            return None

        hf_method = next(
            (method for method in simulation.model_method if isinstance(method, HF)), None
        )
        cc_method = next(
            (method for method in simulation.model_method if isinstance(method, CC)), None
        )
        if hf_method is None or cc_method is None:
            return None

        workflow = SerialWorkflow()
        tasks = [Task(name='HF')]
        if any(
            isinstance(method, OrbitalLocalization) for method in simulation.model_method
        ):
            tasks.append(Task(name='Orbital localization'))
        if any(
            isinstance(method, PerturbationMethod) for method in simulation.model_method
        ):
            tasks.append(Task(name='Local MP2'))

        final_outputs = []
        if simulation.outputs:
            final_outputs.append(Link(name='Outputs', section=simulation.outputs[-1]))
        tasks.append(
            Task(
                name='Local CC' if cc_method.local_correlation else 'CC',
                outputs=final_outputs,
            )
        )
        workflow.tasks = tasks
        self.archive.workflow2 = workflow
        workflow.normalize(self.archive, logger)
        return workflow


class OrcaParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        # Clean up any mapping annotations that might have been left by other parsers
        remove_mapping_annotations(Simulation.m_def)

        reload(orca)

        reader = OutParser()
        reader.filepath = mainfile

        meta = MetainfoParser(data_object=Simulation())
        meta.annotation_key = 'out'
        meta.max_nested_level = 3

        reader.convert(meta)
        archive.data = meta.data_object

        reader.enrich_local_correlation_methods(archive.data)
        reader.enrich_basis_sets(archive.data)
        reader.archive = archive
        archive.workflow2 = reader.build_workflow(archive.data, logger)

        meta.close()
        reader.close()

        # Remove ORCA mapping annotations to avoid interfering with other parsers.
        remove_mapping_annotations(Simulation.m_def)
