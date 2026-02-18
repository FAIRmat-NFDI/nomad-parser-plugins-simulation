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
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

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

        meta.close()
        reader.close()

        # Remove ORCA mapping annotations to avoid interfering with other parsers.
        remove_mapping_annotations(Simulation.m_def)
