from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    pass

import re
import numpy as np
from importlib import reload

from nomad.units import ureg
from nomad.parsing.file_parser import ArchiveWriter, Quantity, TextParser
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path
from nomad_simulation_parsers.parsers.utils.general import remove_mapping_annotations
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulation_parsers.schema_packages import orca
from .text_parser import OutReader
from nomad_simulations.schema_packages import model_system

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
        self.max_nested_level = 5

    def get_program_data(self, src: dict[str, Any]) -> dict[str, Any]:
        return {
            'program_name': 'ORCA',
            'program_version': src.get('program_version'),
        }

    # def get_atoms(self, src: dict[str, Any]):
    #     coords = src.get('single_point', {}).get('cartesian_coordinates', [])
    #     if not coords:
    #         return []

    #     syms, pos = str_to_cartesian_coordinates(coords)
    #     atoms = [{'chemical_symbol': s} for s in syms]
    #     return [{'positions': pos, 'particle_states': atoms}]

    def get_atoms(self, src):
        coords = src.get('single_point', {}).get('cartesian_coordinates', [])
        if not coords:
            return []

        symbols, positions = str_to_cartesian_coordinates(coords)
        n = len(symbols)
        mid = n // 2

        return [{
            'positions'      : positions,
            'particle_states': [{'chemical_symbol': s} for s in symbols],
            'n_particles'    : n,
            'sub_systems'    : [
                {
                    'name'            : 'fragment_A',
                    'type'            : 'molecule / cluster',
                    'particle_indices': list(range(0, mid)),
                    'n_particles'     : mid
                },
                {
                    'name'            : 'fragment_B',
                    'type'            : 'molecule / cluster',
                    'particle_indices': list(range(mid, n)),
                    'n_particles'     : n - mid
                }
            ]
        }]


    
    def get_dft_data(self, source: dict[str, Any]) -> dict[str, Any]:
        """
        Extracts DFT-related data, including XC functionals and SCF settings.
        """
        dft_data = source.get('single_point', {}).get('self_consistent', {}).get('scf_settings', {})
        xc_functionals = []

        # Exchange functional
        if dft_data.get('exchange_functional'):
            xc_functionals.append({
                'libxc_name': dft_data.get('exchange_functional'),
                'name': 'exchange',
                'weight': dft_data.get('scaling_exchange')
            })

        # Correlation functional
        if dft_data.get('correlation_functional'):
            xc_functionals.append({
                'libxc_name': dft_data.get('correlation_functional'),
                'name': 'correlation',
                'weight': dft_data.get('scaling_correlation')
            })
        #print(xc_functionals)
        return {
            'jacobs_ladder': 'metaGGA', # fix here later
            'xc_functionals': xc_functionals,
            'exact_exchange_mixing_factor': dft_data.get('fraction_hf_exchange'),
        }
    
   
class OrcaParser(MatchingParser):
    """
    Minimal NOMAD parser for ORCA.
    """

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
        meta.max_nested_level = 10

        reader.convert(meta)
        archive.data = meta.data_object

        remove_mapping_annotations(orca.general.Simulation.m_def)
        meta.close()
        reader.close()
