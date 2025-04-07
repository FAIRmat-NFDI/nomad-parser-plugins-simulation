import numpy as np
import phonopy
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import (
    AtomicCell,
    AtomsState,
    ModelSystem,
)
from nomad_simulations.schema_packages.outputs import Outputs
from phonopy import Phonopy
from phonopy.units import THzToEv
from structlog.stdlib import BoundLogger

from .calculator import PhononProperties


def write_bandstructure(properties: PhononProperties):
    freqs, bands, bands_labels = properties.get_bandstructure()
    if freqs is None:
        return

    # convert THz to eV
    freqs = freqs * THzToEv

    # convert eV to J
    freqs = (freqs * ureg.eV).to('joules').magnitude

def write_dos(properties: PhononProperties):
    f, dos = properties.get_dos()

    # convert THz to eV to Joules
    f = f * THzToEv
    f = (f * ureg.eV).to('joules').magnitude

def write_thermodynamical_properties(properties: PhononProperties):
    T, fe, _, cv = properties.get_thermodynamical_properties()

    n_atoms = len(properties.phonopy_obj.unitcell)
    n_atoms_supercell = len(properties.phonopy_obj.supercell)

    fe = fe / n_atoms

    # The thermodynamic properties are reported by phonopy for the base
    # system. Since the values in the metainfo are stored per the referenced
    # system, we need to multiple by the size factor between the base system
    # and the supersystem used in the calculations.
    cv = cv * (n_atoms_supercell / n_atoms)

    # convert to SI units
    fe = (fe * ureg.eV).to('joules').magnitude

    cv = (cv * ureg.eV / ureg.K).to('joules/K').magnitude

def create_system(
    cell: np.ndarray,
    symbols: list[str],
    positions: np.ndarray,
    supercell: np.ndarray = None,
) -> ModelSystem:
    sec_system = ModelSystem()
    sec_cell = AtomicCell()
    sec_system.cell.append(sec_cell)

    sec_cell.periodic_boundary_conditions = [True, True, True]
    for symbol in symbols:
        sec_cell.atoms_state.append(AtomsState(chemical_symbol=symbol))

    sec_cell.positions = positions * ureg.angstrom
    sec_cell.lattice_vectors = cell * ureg.angstrom
    sec_cell.supercell_matrix = supercell
    return sec_system


def phonopy_obj_to_archive(
    phonopy_obj: Phonopy,
    archive: EntryArchive = None,
    logger: BoundLogger = None,
    **kwargs,
):
    """
    Run phonopy with an input phonopy object and write the results on a nomad archive.
    """

    logger = logger if logger is not None else get_logger(__name__)
    archive = archive if archive else EntryArchive()

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
    properties = PhononProperties(phonopy_obj, logger, **kwargs)
    properties.get_dos()
    # TODO write outputs
    # vol = np.dot(unit_cell[0], np.cross(unit_cell[1], unit_cell[2]))
    # n_imaginary = np.count_nonzero(properties.frequencies < 0)

    return archive


class PhonopyParser(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger):
        pass
