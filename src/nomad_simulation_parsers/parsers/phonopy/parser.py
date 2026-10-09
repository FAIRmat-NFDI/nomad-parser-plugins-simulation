import os
from typing import Any

import numpy as np
import phonopy
import yaml
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MappingParser, MetainfoParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_system import (
    AlternativeRepresentation,
    AtomsState,
    ModelSystem,
)
from nomad_simulations.schema_packages.variables import Energy2, Temperature
from phonopy.physical_units import get_physical_units
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages.phonopy import (
    PHONOPY_KEY,
    ForceConstants,
    PhononBandStructure,
    PhononDensityOfStates,
    PhonopyMethod,
    PhonopyOutputs,
    PhonopySimulation,
    QPath,
    VibrationalFreeEnergy,
    VibrationalHeatCapacity,
)

from .calculator import PhononProperties

THzToEv = get_physical_units().THzToEv


class PhonopyYamlParser(MappingParser):
    """Read a Phonopy YAML file into the mapping parser's dictionary form."""

    _phonon_properties: PhononProperties | None = None
    _thermodynamic_properties: dict[str, list[dict[str, Any]]] | None = None

    def load_file(self) -> Any:
        if not self.filepath:
            return {}
        phonopy_obj = None
        self._phonon_properties = None
        self._thermodynamic_properties = None
        try:
            mainfile_dir = os.path.dirname(self.filepath)
            force_constants_file = os.path.join(mainfile_dir, 'force_constants.hdf5')
            kwargs = (
                {'force_constants_filename': force_constants_file}
                if os.path.isfile(force_constants_file)
                else {}
            )
            cwd = os.getcwd()
            os.chdir(mainfile_dir)
            try:
                phonopy_obj = phonopy.load(self.filepath, **kwargs)
            finally:
                os.chdir(cwd)
        except Exception:
            return None

        try:
            self._phonon_properties = PhononProperties(phonopy_obj, self.logger)
        except Exception as exc:
            self.logger.warning('Error initializing phonon properties.', exc_info=exc)
        return phonopy_obj

    def to_dict(self) -> dict[str, Any]:
        # Keep the source representation declarative.  The Phonopy object is
        # retained by ``load_file`` for the transform methods below, while
        # the mapping engine receives the YAML dictionary as its source.
        if self.filepath:
            try:
                with open(self.filepath, encoding='utf-8') as handle:
                    return yaml.safe_load(handle) or {}
            except Exception:
                self.logger.error('Error reading phonopy YAML data.')
        return self.data_object if isinstance(self.data_object, dict) else {}

    def get_program(self) -> dict[str, Any]:
        return {'name': 'Phonopy', 'version': phonopy.__version__}

    def get_model_system(self) -> list[dict[str, Any]]:
        obj = self.data_object
        if obj is None:
            return []

        def system(atoms, supercell_matrix=None):
            return {
                'positions': atoms.positions * ureg.angstrom,
                'particle_states': [
                    {'chemical_symbol': symbol} for symbol in atoms.symbols
                ],
                'representations': [
                    {
                        'lattice_vectors': atoms.cell * ureg.angstrom,
                        'periodic_boundary_conditions': [True, True, True],
                        'supercell_matrix': supercell_matrix,
                    }
                ],
            }

        return [system(obj.unitcell), system(obj.supercell, obj.supercell_matrix)]

    def get_model_method(self) -> list[dict[str, Any]]:
        obj = self.data_object
        if obj is None:
            return []
        properties = self._phonon_properties
        if properties is None:
            return []
        displacement = None
        try:
            displacement = np.linalg.norm(obj.displacements[0][1:]) * ureg.angstrom
        except Exception:
            pass
        volume = abs(np.linalg.det(obj.unitcell.cell)) * ureg.angstrom**3
        return [
            {
                'name': 'harmonic lattice dynamics',
                'type': 'finite displacement',
                'displacement': displacement,
                'symmetry_tolerance': obj.symmetry.tolerance * ureg.angstrom,
                'force_calculator': obj.calculator,
                'mesh_density': np.prod(properties.mesh) / volume,
                'with_non_analytic_correction': bool(obj.nac_params),
                'contributions': [],
            }
        ]

    def get_force_constants(self) -> list[dict[str, Any]]:
        obj = self.data_object
        if obj is None:
            return []
        force_constants = obj.force_constants
        if force_constants is None:
            obj.produce_force_constants()
            force_constants = obj.force_constants
        return [{'value': force_constants * ureg.eV / ureg.angstrom**2}]

    def get_phonon_band_structures(
        self,
    ) -> list[dict[str, Any]]:
        self.data_object
        properties = self._phonon_properties
        if properties is None:
            return []
        return [
            {
                'n_bands': segment['frequencies'].shape[-1],
                'value': segment['frequencies'] * ureg.joule,
                'endpoints_labels': [str(label) for label in segment['labels']],
                'q_path': {'points': segment['kpoints']},
            }
            for segment in get_bandstructures(properties)
        ]

    def get_phonon_dos(self) -> list[dict[str, Any]]:
        self.data_object
        properties = self._phonon_properties
        if properties is None:
            return []
        return [
            {
                'value': result['dos'] / ureg.joule,
                'frequencies': {'points': result['frequencies'] * ureg.joule},
            }
            for result in get_dos(properties)
        ]

    def get_vibrational_thermodynamics(
        self,
        key: str | None = None,
    ) -> dict[str, list[dict[str, Any]]] | list[dict[str, Any]]:
        self.data_object
        if self._thermodynamic_properties is not None:
            return (
                self._thermodynamic_properties
                if key is None
                else self._thermodynamic_properties[key]
            )
        properties = self._phonon_properties
        if properties is None:
            self._thermodynamic_properties = {
                'vibrational_free_energies': [],
                'vibrational_heat_capacities': [],
            }
            return (
                self._thermodynamic_properties
                if key is None
                else self._thermodynamic_properties[key]
            )
        thermodynamics = get_thermodynamic_properties(properties)
        temperatures = np.array([entry['temperature'] for entry in thermodynamics])
        temperature_data = {'points': temperatures * ureg.kelvin}
        self._thermodynamic_properties = {
            'vibrational_free_energies': [
                {
                    'value': np.array(
                        [entry['free_energy'] for entry in thermodynamics]
                    )
                    * ureg.joule,
                    'temperatures': temperature_data,
                }
            ],
            'vibrational_heat_capacities': [
                {
                    'value': np.array(
                        [entry['heat_capacity'] for entry in thermodynamics]
                    )
                    * ureg.joule
                    / ureg.kelvin,
                    'temperatures': temperature_data,
                }
            ],
        }
        return (
            self._thermodynamic_properties
            if key is None
            else self._thermodynamic_properties[key]
        )

    def get_imaginary_frequencies(self) -> int:
        self.data_object
        properties = self._phonon_properties
        if properties is None:
            return 0
        if not hasattr(properties, 'frequencies'):
            properties.get_dos()
        return int(np.count_nonzero(properties.frequencies < 0))

    def from_dict(self, dct: dict[str, Any]) -> None:
        # This parser is used as a source parser; it has no mutable target
        # representation to reconstruct.
        self._data = dct


# Compatibility alias for the original misspelling.
PhononyYamlParser = PhonopyYamlParser


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

    # Convert the frequency axis from THz to energy and transform the DOS so
    # that its integral (number of modes) is preserved on the joule axis.
    joules_per_thz = (THzToEv * ureg.eV).to('joules').magnitude
    freq = freq * joules_per_thz
    dos = dos / joules_per_thz
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


def _populate_bandstructures(
    outputs: PhonopyOutputs, properties: PhononProperties
) -> None:
    for segment in get_bandstructures(properties):
        frequencies = segment['frequencies']
        outputs.phonon_band_structures.append(
            PhononBandStructure(
                n_bands=frequencies.shape[-1],
                value=frequencies * ureg.joule,
                endpoints_labels=[str(label) for label in segment['labels']],
                q_path=QPath(points=segment['kpoints']),
            )
        )


def _populate_dos(outputs: PhonopyOutputs, properties: PhononProperties) -> None:
    for dos_result in get_dos(properties):
        outputs.phonon_dos.append(
            PhononDensityOfStates(
                value=dos_result['dos'] / ureg.joule,
                frequencies=Energy2(points=dos_result['frequencies'] * ureg.joule),
            )
        )
    outputs.n_imaginary_frequencies = int(np.count_nonzero(properties.frequencies < 0))


def _populate_thermodynamics(
    outputs: PhonopyOutputs, properties: PhononProperties
) -> None:
    thermodynamics = get_thermodynamic_properties(properties)
    temperatures = np.array([entry['temperature'] for entry in thermodynamics])
    outputs.vibrational_free_energies.append(
        VibrationalFreeEnergy(
            value=np.array([entry['free_energy'] for entry in thermodynamics])
            * ureg.joule,
            temperatures=Temperature(points=temperatures * ureg.kelvin),
        )
    )
    outputs.vibrational_heat_capacities.append(
        VibrationalHeatCapacity(
            value=np.array([entry['heat_capacity'] for entry in thermodynamics])
            * ureg.joule
            / ureg.kelvin,
            temperatures=Temperature(points=temperatures * ureg.kelvin),
        )
    )


def _create_method(
    phonopy_obj: phonopy.Phonopy,
    displacement,
    unit_cell: np.ndarray,
    mesh: list[int],
) -> PhonopyMethod:
    volume = abs(np.linalg.det(unit_cell)) * ureg.angstrom**3
    return PhonopyMethod(
        name='harmonic lattice dynamics',
        type='finite displacement',
        displacement=displacement,
        symmetry_tolerance=phonopy_obj.symmetry.tolerance * ureg.angstrom,
        force_calculator=phonopy_obj.calculator,
        mesh_density=np.prod(mesh) / volume,
        with_non_analytic_correction=bool(phonopy_obj.nac_params),
    )


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

    unit_cell = phonopy_obj.unitcell.cell
    unit_pos = phonopy_obj.unitcell.positions
    unit_sym = phonopy_obj.unitcell.symbols

    super_cell = phonopy_obj.supercell.cell
    super_pos = phonopy_obj.supercell.positions
    super_sym = phonopy_obj.supercell.symbols

    try:
        displacement = np.linalg.norm(phonopy_obj.displacements[0][1:])
        displacement = displacement * ureg.angstrom
    except Exception:
        displacement = None

    supercell_matrix = phonopy_obj.supercell_matrix
    data = Simulation()
    archive.data = data

    data.program = Program(name='Phonopy', version=phonopy.__version__)

    sec_system_unit = create_system(unit_cell, unit_sym, unit_pos)
    data.model_system.append(sec_system_unit)

    sec_system = create_system(super_cell, super_sym, super_pos, supercell_matrix)
    data.model_system.append(sec_system)

    try:
        force_constants = phonopy_obj.force_constants
        if force_constants is None:
            phonopy_obj.produce_force_constants()
            force_constants = phonopy_obj.force_constants
        force_constants = force_constants * ureg.eV / ureg.angstrom**2
    except Exception:
        logger.error('Error producing force constants.')
        return

    properties = PhononProperties(phonopy_obj, logger, **kwargs)
    method = _create_method(phonopy_obj, displacement, unit_cell, properties.mesh)
    data.model_method.append(method)

    sec_outputs = PhonopyOutputs()
    sec_outputs.model_system_ref = sec_system
    sec_outputs.model_method_ref = method
    sec_outputs.force_constants.append(ForceConstants(value=force_constants))
    data.outputs.append(sec_outputs)

    # run Phonopy
    try:
        _populate_bandstructures(sec_outputs, properties)
    except Exception as exc:
        logger.warning('Error calculating phonon band structure.', exc_info=exc)

    try:
        _populate_dos(sec_outputs, properties)
    except Exception as exc:
        logger.warning('Error calculating phonon density of states.', exc_info=exc)

    try:
        _populate_thermodynamics(sec_outputs, properties)
    except Exception as exc:
        logger.warning('Error calculating phonon thermodynamics.', exc_info=exc)

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
                positions=atoms.positions * ureg.angstrom,
                cell=atoms.cell * ureg.angstrom,
                particle_states=[dict(chemical_symbol=sym) for sym in atoms.symbols],
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
        mainfile = os.path.abspath(self.mainfile)
        mainfile_dir = os.path.dirname(mainfile)
        cwd = os.getcwd()
        os.chdir(mainfile_dir)
        try:
            # ``phonopy.yaml`` commonly stores the structure and references a
            # neighbouring HDF5 force-constant file.  ``phonopy.load`` does
            # not discover that sidecar automatically, so make the parser's
            # public entry point behave like the usual phonopy directory
            # layout.
            force_constants_file = os.path.join(mainfile_dir, 'force_constants.hdf5')
            if os.path.isfile(force_constants_file):
                phonopy_obj = phonopy.load(
                    mainfile,
                    force_constants_filename=force_constants_file,
                )
            else:
                phonopy_obj = phonopy.load(mainfile)
        except Exception:
            self.logger.error('Error loading phonopy file.')
            phonopy_obj = None
        finally:
            os.chdir(cwd)

        if phonopy_obj is None:
            return

        phonopy_obj_to_archive(phonopy_obj, self.archive, self.logger)


class PhonopyMappingArchiveWriter(PhonopyArchiveWriter):
    """Archive writer using annotation-driven mapping into the simulation schema."""

    def write_to_archive(self):
        mainfile = os.path.abspath(self.mainfile)
        source_parser = PhonopyYamlParser(filepath=mainfile, logger=self.logger)
        if not source_parser.data:
            self.logger.error('Error parsing phonopy YAML file.')
            return

        parser = MetainfoParser(data_object=PhonopySimulation(), logger=self.logger)
        parser.annotation_key = PHONOPY_KEY
        source_parser.convert(parser)
        self.archive.data = parser.data_object

        # References are intentionally linked after conversion: the generic
        # mapper cannot resolve object references between freshly-created list
        # elements from a plain source dictionary.
        if self.archive.data.outputs:
            self.archive.data.outputs[
                0
            ].model_system_ref = self.archive.data.model_system[-1]
            self.archive.data.outputs[
                0
            ].model_method_ref = self.archive.data.model_method[0]


class PhonopyParser(MatchingParser):
    archive_writer = PhonopyMappingArchiveWriter()

    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger):
        self.archive_writer.write(mainfile, archive, logger)
