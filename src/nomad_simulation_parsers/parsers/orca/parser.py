from typing import TYPE_CHECKING, Any

import numpy as np
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad_simulations.schema_packages.atoms_state import AtomsState
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem
from nomad_simulations.schema_packages.outputs import Outputs
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
)

from .text_parser import OutReader

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import EntryArchive
    from structlog.stdlib import BoundLogger

CARTESIAN_COORDINATE_LENGTH = 3
CARTESIAN_COORDINATE_STRIDE = 4
ORBITAL_ENERGY_N_COLUMNS = 4
ORBITAL_ENERGY_N_DIMENSIONS = 2


def _as_dict(value: Any) -> dict[str, Any]:
    return value._results if hasattr(value, '_results') else value or {}


def _first(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def _cartesian_coordinates_to_system(
    coordinates: list[Any] | None,
    scf_settings: dict[str, Any],
) -> ModelSystem | None:
    if not coordinates:
        return None

    symbols = []
    positions = []
    for index in range(0, len(coordinates), CARTESIAN_COORDINATE_STRIDE):
        symbol = coordinates[index]
        coordinate = coordinates[index + 1 : index + CARTESIAN_COORDINATE_STRIDE]
        if (
            not isinstance(symbol, str)
            or len(coordinate) != CARTESIAN_COORDINATE_LENGTH
        ):
            continue
        symbols.append(symbol)
        positions.append(coordinate)

    if not symbols:
        return None

    system = ModelSystem(
        is_representative=True,
        positions=np.asarray(positions, dtype=np.float64) * ureg.angstrom,
        particle_states=[AtomsState(chemical_symbol=symbol) for symbol in symbols],
    )

    total_charge = _first(scf_settings.get('total_charge'))
    if total_charge is not None:
        system.total_charge = int(total_charge)

    multiplicity = _first(scf_settings.get('multiplicity'))
    if multiplicity is not None:
        system.total_spin = int(multiplicity) - 1

    return system


def _build_molecular_orbitals(
    orbital_energies: list[np.ndarray] | np.ndarray | None,
    n_ao: int | None,
) -> MolecularOrbitals | None:
    if orbital_energies is None:
        return None

    table = _first(orbital_energies)
    if table is None:
        return None

    table = np.asarray(table, dtype=np.float64)
    if (
        table.ndim != ORBITAL_ENERGY_N_DIMENSIONS
        or table.shape[1] < ORBITAL_ENERGY_N_COLUMNS
    ):
        return None

    return MolecularOrbitals(
        n_mo=table.shape[0],
        n_ao=n_ao,
        mo_occupations=table[:, 1],
        mo_energies=table[:, 3] * ureg.electron_volt,
        mo_type='canonical',
    )


class OrcaParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        reader = OutReader()
        reader.mainfile = mainfile
        reader.parse()

        results = reader.results or {}
        single_point = _as_dict(results.get('single_point'))
        self_consistent = _as_dict(single_point.get('self_consistent'))
        scf_settings = _as_dict(self_consistent.get('scf_settings'))
        basis_set_total = _as_dict(results.get('basis_set_total'))

        simulation = Simulation(
            program=Program(name='ORCA', version=results.get('program_version'))
        )

        system = _cartesian_coordinates_to_system(
            single_point.get('cartesian_coordinates'), scf_settings
        )
        if system is not None:
            simulation.model_system = [system]

        molecular_orbitals = _build_molecular_orbitals(
            self_consistent.get('orbital_energies'),
            n_ao=_first(basis_set_total.get('main_basis_set')),
        )
        if molecular_orbitals is not None:
            outputs = Outputs(electronic_eigenvalues=[molecular_orbitals])
            if system is not None:
                outputs.model_system_ref = system
            simulation.outputs = [outputs]

        archive.data = simulation
