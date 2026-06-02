from importlib import reload
from typing import TYPE_CHECKING, Any

import numpy as np
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
)

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

    def get_molecular_orbitals(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        single_point = self._as_dict(src.get('single_point', {}))
        self_consistent = self._as_dict(single_point.get('self_consistent', {}))
        basis_set_total = self._as_dict(src.get('basis_set_total', {}))

        orbital_energies = self._scalar(self_consistent.get('orbital_energies'))
        if orbital_energies is None:
            return []

        table = np.asarray(orbital_energies, dtype=np.float64)
        if (
            table.ndim != ORBITAL_ENERGY_N_DIMENSIONS
            or table.shape[1] < ORBITAL_ENERGY_N_COLUMNS
        ):
            return []

        n_ao = self._scalar(basis_set_total.get('main_basis_set'))
        return [
            {
                'm_def': MolecularOrbitals.m_def.qualified_name(),
                'n_mo': int(table.shape[0]),
                'n_ao': int(n_ao) if n_ao is not None else int(table.shape[0]),
                'mo_occupations': table[:, ORBITAL_OCCUPATION_COLUMN],
                'mo_energies': table[:, ORBITAL_ENERGY_EV_COLUMN] * ureg.electron_volt,
                'mo_type': 'canonical',
            }
        ]

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
        metainfo_parser = MetainfoParser(data_object=Simulation())
        metainfo_parser.annotation_key = orca.OUT_KEY
        metainfo_parser.max_nested_level = 3

        try:
            reader.convert(metainfo_parser)
            archive.data = metainfo_parser.data_object
        finally:
            metainfo_parser.close()
            reader.close()
