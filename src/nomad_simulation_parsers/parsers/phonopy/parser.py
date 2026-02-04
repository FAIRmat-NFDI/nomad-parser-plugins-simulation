import os
from typing import Any

import numpy as np
import phonopy
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import (
    AlternativeRepresentation,
    AtomsState,
    ModelSystem,
)
from nomad_simulations.schema_packages.outputs import Outputs
from phonopy.units import THzToEv
from structlog.stdlib import BoundLogger

from .calculator import PhononProperties


def get_bandstructures(properties: PhononProperties) -> list[dict[str, Any]]:
    freqs, bands, bands_labels = properties.get_bandstructure()
    if freqs is None:
        return []

    # convert THz to eV
    freqs = freqs * THzToEv

    # convert eV to J
    freqs = (freqs * ureg.eV).to('joules').magnitude
    return [
        dict(frequencies=freq, kpoints=bands[n], labels=bands_labels[n])
        for n, freq in enumerate(freqs)
    ]


def get_dos(properties: PhononProperties) -> list[dict[str, Any]]:
    freq, dos = properties.get_dos()

    # convert THz to eV to Joules
    freq = freq * THzToEv
    freq = (freq * ureg.eV).to('joules').magnitude
    return [dict(frequencies=freq, dos=dos)]


def get_thermodynamic_properties(properties: PhononProperties) -> list[dict[str, Any]]:
    temperatures, free_energies, _, heat_capacties = (
        properties.get_thermodynamical_properties()
    )
    n_atoms = len(properties.phonopy_obj.unitcell)
    n_atoms_supercell = len(properties.phonopy_obj.supercell)

    free_energies = free_energies / n_atoms

    # The thermodynamic properties are reported by phonopy for the base
    # system. Since the values in the metainfo are stored per the referenced
    # system, we need to multiple by the size factor between the base system
    # and the supersystem used in the calculations.
    heat_capacties = heat_capacties * (n_atoms_supercell / n_atoms)

    # convert to SI units
    free_energies = (free_energies * ureg.eV).to('joules').magnitude

    heat_capacties = (heat_capacties * ureg.eV / ureg.K).to('joules/K').magnitude
    return [
        dict(
            temperature=temperature,
            free_energy=free_energies[n],
            heat_capacity=heat_capacties[n],
        )
        for n, temperature in enumerate(temperatures)
    ]


def create_system(
    cell: np.ndarray,
    symbols: list[str],
    positions: np.ndarray,
    supercell: np.ndarray = None,
) -> ModelSystem:
    sec_system = ModelSystem()
    sec_representation = AlternativeRepresentation()
    sec_system.representations.append(sec_representation)

    sec_representation.periodic_boundary_conditions = [True, True, True]
    for symbol in symbols:
        sec_system.particle_states.append(AtomsState(chemical_symbol=symbol))

    sec_system.positions = positions * ureg.angstrom
    sec_representation.lattice_vectors = cell * ureg.angstrom
    sec_representation.supercell_matrix = supercell
    return sec_system


def phonopy_obj_to_archive(
    phonopy_obj: phonopy.Phonopy,
    archive: EntryArchive = None,
    logger: BoundLogger = None,
    **kwargs,
):
    """
    Run phonopy with an input phonopy object and write the results on a nomad archive.
    """

    logger = logger if logger is not None else get_logger(__name__)
    archive = archive if archive is not None else EntryArchive()

    unit_cell = phonopy_obj.unitcell.get_cell()
    unit_pos = phonopy_obj.unitcell.get_positions()
    unit_sym = phonopy_obj.unitcell.get_chemical_symbols()

    super_cell = phonopy_obj.supercell.get_cell()
    super_pos = phonopy_obj.supercell.get_positions()
    super_sym = phonopy_obj.supercell.get_chemical_symbols()

    try:
        displacement = np.linalg.norm(phonopy_obj.displacements[0][1:])
        displacement = displacement * ureg.angstrom
    except Exception:
        displacement = None

    supercell_matrix = phonopy_obj.supercell_matrix
    # sym_tol = phonopy_obj.symmetry.tolerance

    data = Simulation()
    archive.data = data

    data.program = Program(name='Phonopy', version=phonopy.__version__)

    sec_system_unit = create_system(unit_cell, unit_sym, unit_pos)
    data.model_system.append(sec_system_unit)

    sec_system = create_system(super_cell, super_sym, super_pos, supercell_matrix)
    data.model_system.append(sec_system)

    try:
        force_constants = phonopy_obj.get_force_constants()
        force_constants = (
            (force_constants * ureg.eV / ureg.angstrom**2).to('J/(m**2)').magnitude
        )
    except Exception:
        logger.error('Error producing force constants.')
        return

    sec_outputs = Outputs()
    sec_outputs.model_system_ref = sec_system

    # run Phonopy
    # properties = PhononProperties(phonopy_obj, logger, **kwargs)
    # TODO write outputs

    return archive


def phonopy_obj_to_dict(
    phonopy_obj: phonopy.Phonopy, logger: BoundLogger = None, **kwargs
) -> dict[str, Any]:
    logger = logger if logger is not None else get_logger(__name__)

    results = dict(program=dict(name='Phonopy', version=phonopy.__version__))

    system = results.setdefault('model_system', [])
    for atoms in [phonopy_obj.unitcell, phonopy_obj.supercell]:
        system.append(
            dict(
                positions=atoms.get_positions() * ureg.angstrom,
                cell=atoms.get_cell() * ureg.angstrom,
                particle_states=[
                    dict(chemical_symbol=sym) for sym in atoms.get_chemical_symbols()
                ],
            )
        )
    if system:
        system[-1]['supercell_matrix'] = phonopy_obj.supercell_matrix

    # run Phonopy
    properties = PhononProperties(phonopy_obj, logger, **kwargs)
    results['outputs'] = dict(
        dos=get_dos(properties),
        bandstructures=get_bandstructures(properties),
        thermodynamics=get_thermodynamic_properties(properties),
    )

    return results


class PhonopyArchiveWriter(ArchiveWriter):
    def write_to_archive(self):
        cwd = os.getcwd()
        os.chdir(os.path.dirname(self.mainfile))
        try:
            phonopy_obj = phonopy.load(self.mainfile)
        except Exception:
            self.logger.error('Error loading phonopy file.')
            phonopy_obj = None
        finally:
            os.chdir(cwd)

        if phonopy_obj is None:
            return

        phonopy_obj_to_archive(phonopy_obj, self.archive, self.logger)


class PhonopyParser(MatchingParser):
    archive_writer = PhonopyArchiveWriter()

    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger):
        self.archive_writer.write(mainfile, archive, logger)
