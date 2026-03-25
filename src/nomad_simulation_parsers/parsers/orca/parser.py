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
        Minimal DFT getter

        """
        scf_settings = (
            src.get('single_point', {})
            .get('self_consistent', {})
            .get('scf_settings', {})
        )
        if not scf_settings:
            return {}

        hf_frac = scf_settings.get('fraction_hf_exchange')

        return {'xc': {'global_exact_exchange': hf_frac}}

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

        meta = MetainfoParser(data_object=Simulation())
        meta.annotation_key = 'out'
        # meta.max_nested_level = 1

        reader.convert(meta)
        archive.data = meta.data_object

        meta.close()
        reader.close()
