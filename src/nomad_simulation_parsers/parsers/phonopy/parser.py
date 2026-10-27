import os
from typing import Any

import numpy as np
import phonopy
import yaml
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MappingParser, MetainfoParser
from phonopy.interface.calculator import (
    get_default_physical_units,
    get_force_constant_conversion_factor,
)
from phonopy.physical_units import get_physical_units
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages.phonopy import (
    PHONOPY_KEY,
    PhonopySimulation,
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

    def get_model_ref(self, key: str = 'system') -> str | None:
        obj = self.data_object
        if obj is None:
            return
        if key == 'system' and obj.unitcell:
            return '/data/model_system/1'
        if key == 'method' and self._phonon_properties is not None:
            return '/data/model_method/0'
        return None

    def get_model_system(self) -> list[dict[str, Any]]:
        obj = self.data_object
        if obj is None:
            return []

        return [
            {
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
            for atoms, supercell_matrix in (
                (obj.unitcell, None),
                (obj.supercell, obj.supercell_matrix),
            )
        ]

    def get_model_method(self) -> list[dict[str, Any]]:
        obj = self.data_object
        if obj is None:
            return []
        properties = self._phonon_properties
        if properties is None:
            return []
        displacement = None
        try:
            first_displacement = np.asarray(obj.displacements[0])
            if first_displacement.ndim == 1:
                first_displacement = first_displacement[1:]
            displacement = np.linalg.norm(first_displacement) * ureg.angstrom
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
            try:
                obj.produce_force_constants()
            except Exception as exc:
                self.logger.warning('Error producing force constants.', exc_info=exc)
                return []
            force_constants = obj.force_constants
        if force_constants is None:
            return []
        native_unit = get_default_physical_units(obj.calculator)['force_constants_unit']
        conversion_factor = get_force_constant_conversion_factor(native_unit, 'vasp')
        return [
            {'value': force_constants * conversion_factor * ureg.eV / ureg.angstrom**2}
        ]

    def get_phonon_band_structures(
        self,
    ) -> list[dict[str, Any]]:
        self.data_object
        properties = self._phonon_properties
        if properties is None:
            return []
        frequencies, bands, labels = properties.get_bandstructure()
        if frequencies is None:
            return []
        frequencies = (frequencies * THzToEv * ureg.eV).to('joule').magnitude
        return [
            {
                'n_bands': segment.shape[-1],
                'value': segment * ureg.joule,
                'endpoints_labels': [str(label) for label in labels[index]],
                'q_path': {'points': bands[index]},
            }
            for index, segment in enumerate(frequencies)
        ]

    def get_phonon_dos(self) -> list[dict[str, Any]]:
        self.data_object
        properties = self._phonon_properties
        if properties is None:
            return []
        frequencies, dos = properties.get_dos()
        joules_per_thz = (THzToEv * ureg.eV).to('joule').magnitude
        return [
            {
                'value': dos / joules_per_thz / ureg.joule,
                'frequencies': {'points': frequencies * joules_per_thz * ureg.joule},
            }
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
        temperatures, free_energies, _, heat_capacities = (
            properties.get_thermodynamical_properties()
        )
        n_atoms = len(properties.phonopy_obj.unitcell)
        n_atoms_supercell = len(properties.phonopy_obj.supercell)
        free_energies = (free_energies / n_atoms * ureg.eV).to('joule').magnitude
        heat_capacities = (
            (heat_capacities * (n_atoms_supercell / n_atoms) * ureg.eV / ureg.kelvin)
            .to('joule / kelvin')
            .magnitude
        )
        temperature_data = {'points': temperatures * ureg.kelvin}
        self._thermodynamic_properties = {
            'vibrational_free_energies': [
                {
                    'value': free_energies * ureg.joule,
                    'temperatures': temperature_data,
                }
            ],
            'vibrational_heat_capacities': [
                {
                    'value': heat_capacities * ureg.joule / ureg.kelvin,
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


class PhonopyArchiveWriter(ArchiveWriter):
    """Archive writer using annotation-driven mapping into the simulation schema."""

    def write(
        self,
        source: str | os.PathLike[str] | phonopy.Phonopy,
        archive: EntryArchive,
        logger: BoundLogger | None = None,
        child_archives=None,
    ) -> None:
        """Write a Phonopy YAML file or an in-memory Phonopy calculation."""
        super().write(source, archive, logger, child_archives)

    def write_to_archive(self):
        if isinstance(self.mainfile, phonopy.Phonopy):
            source_parser = PhonopyYamlParser(
                data_object=self.mainfile, logger=self.logger
            )
            try:
                source_parser._phonon_properties = PhononProperties(
                    self.mainfile, self.logger
                )
            except Exception as exc:
                self.logger.warning(
                    'Error initializing phonon properties.', exc_info=exc
                )
        else:
            mainfile = os.path.abspath(self.mainfile)
            source_parser = PhonopyYamlParser(filepath=mainfile, logger=self.logger)
            if not source_parser.data:
                self.logger.error('Error parsing phonopy YAML file.')
                return

        parser = MetainfoParser(data_object=PhonopySimulation(), logger=self.logger)
        parser.annotation_key = PHONOPY_KEY
        source_parser.convert(parser)
        self.archive.data = parser.data_object

        # close parsers
        source_parser.close()
        parser.close()


class PhonopyParser(MatchingParser):
    archive_writer = PhonopyArchiveWriter()

    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger):
        self.archive_writer.write(mainfile, archive, logger)
