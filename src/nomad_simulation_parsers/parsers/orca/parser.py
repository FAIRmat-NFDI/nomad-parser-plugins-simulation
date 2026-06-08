import re
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
        finally:
            metainfo_parser.close()
            reader.close()
