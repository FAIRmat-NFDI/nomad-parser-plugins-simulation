from collections.abc import Sequence
from importlib import reload
from typing import Any

import numpy as np
import pint
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser.mapping_parser import HDF5Parser, MetainfoParser, Path
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger

# from nomad_simulations.schema_packages.workflow import molecular_dynamics
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser
from nomad_simulation_parsers.schema_packages import h5md
from nomad_simulation_parsers.schema_packages.h5md import MolecularDynamics, Simulation
from nomad_simulation_parsers.schema_packages.utils import remove_mapping_annotations

LOGGER = get_logger(__name__)


class H5MDMetainfoParser(MetainfoParser):
    # TODO: temporary fix until structlog propagation lands everywhere
    @property
    def logger(self):
        return LOGGER


class H5MDH5Parser(HDF5Parser):
    trajectory_steps: list[int] = []

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_value(self, name: str, dct: dict[str, Any]) -> Any:
        value = dct.get(name, {})
        if not isinstance(value, dict):
            return value
        value = value.get(self.value_key)
        if value is None:
            return
        unit = dct.get(name, {}).get(f'{self.attribute_prefix}unit')
        if unit:
            value = value * ureg(unit)
        factor = dct.get(name, {}).get(f'{self.attribute_prefix}unit_factor')
        if factor:
            value = value * factor
        return value

    def get_source(self, parent: dict[str, Any], path: str) -> Any:
        if path is None:
            return {}

        path_segments = path.split('.', 1)
        source = parent.get(path_segments[0], {})

        if len(path_segments) == 1:
            return source
        return self.get_source(source, path_segments[1])

    def map_value(
        self, source: dict[str, Any], key: str = None, enum_spec: str = None
    ) -> Any:
        if key is None:
            return None

        value = self.get_value(key, source)

        enum_map = {'upper': str.upper, 'lower': str.lower}

        if enum_spec in enum_map and isinstance(value, str):
            return enum_map[enum_spec](value)

        return value

    def get_sub_systems(self, source: dict[str, Any], **kwargs) -> list[dict[str, Any]]:
        step = source.get('step', None)
        if step is not None:
            if step != 0:  # TODO extend to time-dependent bond lists and topologies
                return []
            source = self.get_source(self.data, kwargs['path'])

        particles_group = source.get('particles_group', {})

        source = list(group for group in particles_group.values())

        return source

    def get_traj_data(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        positions = self.get_source(source, 'position')
        steps = self.get_value('step', positions)
        times = self.get_value('time', positions)
        positions = self.get_value('value', positions)
        velocities = self.get_source(source, 'velocity')
        velocities = self.get_value('value', velocities)

        if not (len(steps) == len(times) == len(positions)):
            self.logger.error('Trajectory data length mismatch')
            return []

        if velocities is not None and len(positions) != len(velocities):
            self.logger.error('Velocity length mismatch')
            return []
        # TODO add len assertion for other system traj properties
        # TODO or generalize to store properly

        traj_data = [
            {
                'step': step,
                'time': times[n],
                'n_particles': len(positions[n]) if len(positions) != 0 else [],
                'positions': positions[n] if len(positions) != 0 else [],
                'velocities': velocities[n] if len(velocities) != 0 else [],
            }
            for n, step in enumerate(steps)
            if step in self.trajectory_steps
        ]

        return traj_data

    def get_step_data(self, data: dict[str, Any], step: int) -> dict[str, Any]:
        step_data = {}
        value = self.get_value('value', data)
        steps = self.get_value('step', data)
        if value is None or steps is None:
            return step_data
        index = steps.index(step)
        step_data['value'] = value[index]
        times = self.get_value('time', data)
        step_data['time'] = times[index]
        return step_data

    def get_cell_data(self, source: dict[str, Any]) -> dict[str, Any]:
        if not source.get('step'):
            return {}
        particles = self.data.get('particles', {}).get('all')
        if particles is None:
            return {}

        source_paths = [
            ('lattice_vectors', 'box.edges'),
        ]
        system_data = {}
        for name, path in source_paths:
            data = self.get_source(particles, path)
            if not data:
                continue
            step_data = self.get_step_data(data, source.get('step'))
            if not step_data:
                continue
            system_data.setdefault('time', step_data['time'])
            if step_data['time'] == system_data['time']:
                system_data[name] = step_data['value']
        box = particles.get('box', {})
        system_data['boundary'] = box.get(f'{self.attribute_prefix}boundary')

        return system_data

    def to_species_labels(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        if source.get('step') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])

        return [{'chemical_symbol': s, 'label': s} for s in source_data]

    def get_top_system_quantity(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        if source.get('step') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])
        if kwargs.get('func'):
            try:
                source_data = kwargs.get('func')(source_data)
            except Exception:
                self.logger.warning(
                    'Error applying function to get top level system quantity, '
                    'will not be populated.'
                )
                return []

        return source_data

    def get_output_steps(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        output_steps = {}

        def get_steps(dct: dict[str, Any]) -> dict[str, Any]:
            steps = self.get_value('step', dct)
            if steps is None:
                return {}
            times = self.get_value('time', dct)
            if len(steps) != len(times):
                self.logger.error(
                    'Inconsistent step-time combinations in observable data.'
                )

            return {step: times[n] for n, step in enumerate(steps)}

        def get_observable_steps(
            source: dict[str, Any], output_steps: dict[str, Any]
        ) -> list[dict[str, Any]]:
            for __, val in source.items():
                observable_type = val.get('@type')
                if not observable_type:
                    get_observable_steps(val, output_steps)
                elif observable_type == 'configurational':
                    steps = get_steps(val)
                    if not all(
                        [
                            time == output_steps[step]
                            for step, time in steps.items()
                            if step in output_steps.keys()
                        ]
                    ):
                        self.logger.error(
                            'Inconsistent step-time combinations in observable data.'
                        )
                    output_steps.update(steps)

            return [{'step': step, 'time': time} for step, time in output_steps.items()]

        return get_observable_steps(source, output_steps)

    def get_contributions(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        if source.get('step') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])
        include = kwargs.get('include')
        exclude = kwargs.get('exclude')
        contributions = []
        for key, val in source_data.items():
            # Skip if key is specifically excluded
            if exclude and key in exclude:
                continue
            # Skip if include list is provided but key is not in it
            if include and key not in include:
                continue
            step_data = self.get_step_data(val, source['step'])
            contributions.append({'name': key, **step_data})
        return contributions

    def get_output_data(
        self, source: dict[str, Any], **kwargs
    ) -> pint.Quantity | list[dict[str, Any]] | None:
        """Router function that delegates to appropriate output extraction method."""
        if source.get('value') is not None:
            return source['value']

        # Route based on whether this is step-based (configurational) or stepless data
        if source.get('step') is not None:
            # Step-based: configurational outputs (energies, forces, temperatures)
            return self.get_configurational_output(source, **kwargs)
        else:
            # Stepless: ensemble_average or correlation_function outputs (RDFs, MSDs)
            return self.get_ensemble_output(source, **kwargs)

    # TODO reassess the necessity of this function
    def _get_output_data_list_from_source(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        """Extract multiple objects directly from source (e.g., RDFs, MSDs)."""
        try:
            results = []

            # Iterate through each item in the source (e.g., MOL1-MOL1, MOL1-MOL2, etc.)
            for item_name, item_data in source.items():
                # Create entry with label and extract specific values using get_value
                entry = {'label': item_name}

                # Extract common RDF/MSD quantities
                for quantity in ['bins', 'value']:
                    value = self.get_value(quantity, item_data)
                    if value is not None:
                        entry[quantity] = value

                results.append(entry)

            return results

        except Exception as e:
            self.logger.warning(f'Could not extract output data list from source: {e}')
            return []

    def get_custom_outputs(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        """Extract custom configurational outputs (step-based data)."""

        # Only process step-based (configurational) data
        if source.get('step') is None:
            return []

        source_data = self.get_source(self.data, kwargs['path'])
        include = kwargs.get('include')
        exclude = kwargs.get('exclude', [])

        custom_outputs = []
        for key, val in source_data.items():
            # Skip standard observables
            if exclude and key in exclude:
                continue

            # Only include step-based data (configurational custom outputs)
            if val.get('step') is not None:
                if include is None or key in include:
                    step_data = self.get_step_data(val, source['step'])
                    if step_data.get('value') is not None:
                        if isinstance(step_data['value'], pint.Quantity):
                            step_data['unit'] = str(step_data['value'].units)
                            step_data['value'] = step_data['value'].magnitude
                    custom_outputs.append({'name': key, **step_data})

        return custom_outputs

    def get_custom_ensemble_outputs(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        """
        Get custom ensemble observables following the same pattern as
        get_ensemble_output. Filters for observables with type 'ensemble_average'
        or 'correlation_function'. Returns list of dicts for mapper auto-generation
        of subsection members.
        """
        # Use EXACT same logic as get_ensemble_output but with type filtering
        if source.get('value') is not None:
            return source['value']

        # For stepless data, convert to list format if it contains multiple items
        if not source or not all(isinstance(v, dict) for v in source.values()):
            return []

        results = []
        exclude_list = kwargs.get('exclude', [])
        type_filter = kwargs.get(
            'type_filter', ['ensemble_average', 'correlation_function']
        )

        # Convert source = {label_1: dict_1, ...} to results = [dict_1*, ...]
        # where dict_x* = dict_x + {label: label_x}
        for label, data_dict in source.items():
            # Skip standard properties listed in exclude list
            if label in exclude_list:
                continue

            # First check for nested structures (like free_energies with
            # individual states)
            nested_observables_found = False
            if isinstance(data_dict, dict):
                for nested_label, nested_data in data_dict.items():
                    if isinstance(nested_data, dict):
                        nested_type = self._get_observable_type(
                            f'{label}/{nested_label}', nested_data
                        )
                        if nested_type in type_filter:
                            nested_observables_found = True
                            result_dict = {'label': f'{label}_{nested_label}'}
                            # Process nested data with consistent logic
                            self._process_observable_data(
                                result_dict, nested_data, nested_type, type_filter
                            )
                            results.append(result_dict)

            # If no nested observables found, process this level directly
            if not nested_observables_found:
                data_type = self._get_observable_type(label, data_dict)
                if data_type in type_filter:
                    result_dict = {'label': label}
                    # Process data with consistent logic
                    self._process_observable_data(
                        result_dict, data_dict, data_type, type_filter
                    )
                    results.append(result_dict)
        return results

    def _process_observable_data(
        self, result_dict: dict, data_dict: dict, data_type: str, type_filter: list
    ) -> None:
        """Helper method to consistently process observable data."""
        for key in data_dict.keys():
            # Skip @type, None values, and 'times' outside correlation functions
            if (
                key == '@type'
                or (value := self.get_value(key, data_dict)) is None
                or (key == 'times' and data_type != 'correlation_function')
            ):
                continue
            self._add_value_to_result(result_dict, key, value)

    def _add_value_to_result(self, result_dict: dict, key: str, value) -> None:
        """Helper method to add a value to result dict with proper unit handling."""
        # Handle 'times' field as normal Quantity
        if key == 'times':
            result_dict[key] = value
        # Handle 'value' field specially - always use schema-mapped field names
        elif key == 'value':
            if hasattr(value, 'magnitude') and hasattr(value, 'units'):
                # Quantity with units
                magnitude = value.magnitude
                # Convert scalars to 1-element arrays for schema consistency
                if not isinstance(magnitude, np.ndarray | Sequence):
                    magnitude = [magnitude]
                result_dict['value_magnitude'] = magnitude
                result_dict['value_unit'] = str(value.units)
            else:
                # Raw value (list, float, etc.) - still use mapped field name
                # Convert scalars to 1-element arrays for schema consistency
                if not isinstance(value, np.ndarray | Sequence):
                    value = [value]
                result_dict['value_magnitude'] = value
        # Handle other Quantities by separating magnitude and unit
        elif hasattr(value, 'magnitude') and hasattr(value, 'units'):
            result_dict[key] = value.magnitude
            result_dict[f'{key}_unit'] = str(value.units)
        else:
            result_dict[key] = value

    def _get_observable_type(self, label: str, data_dict: dict) -> str | None:
        """Get the type of an observable from H5MD file or parsed data."""
        # Try to get type from raw H5MD file attributes
        if hasattr(self, 'h5_parser') and hasattr(self.h5_parser, 'h5_archive'):
            try:
                # Access h5_archive directly (already opened by h5_parser)
                # Will be closed via h5_parser.close()
                h5_file = self.h5_parser.h5_archive
                if h5_file and 'observables' in h5_file:
                    obs_path = f'observables/{label}'
                    if obs_path in h5_file:
                        obs_group = h5_file[obs_path]
                        if hasattr(obs_group, 'attrs') and 'type' in obs_group.attrs:
                            return obs_group.attrs['type']
            except Exception:
                pass
            finally:
                self.h5_parser.h5_archive.close()

        # Fallback to parsed data structure
        return data_dict.get('@type') or data_dict.get('attrs', {}).get('type')

    def get_configurational_output(
        self, source: dict[str, Any], **kwargs
    ) -> pint.Quantity | None:
        """Extract configurational outputs (energies, forces, temperatures)."""
        if source.get('value') is not None:
            return source['value']

        # This function is only for step-based (configurational) outputs
        if source.get('step') is None:
            return None

        source_data = self.get_source(self.data, kwargs['path'])
        if not source_data:
            self.logger.warning(
                'No data found at specified path for observable extraction.'
            )
            return None

        # For configurational outputs, we can optionally validate @type if present
        actual_type = source_data.get('@type')
        if actual_type is not None and actual_type != 'configurational':
            self.logger.warning(
                'Expected configurational output but found different type. '
                'Skipping this observable.'
            )
            return None
        elif actual_type is None:
            self.logger.debug(
                'No type defined in source data, proceeding with extraction.'
            )

        step_data = self.get_step_data(source_data, source['step'])
        return step_data.get('value')

    def get_ensemble_output(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]] | None:
        """Extract ensemble/correlation outputs (RDFs, MSDs) - stepless data."""
        if source.get('value') is not None:
            return source['value']

        # For stepless data, convert to list format if it contains multiple items
        if source and all(isinstance(v, dict) for v in source.values()):
            results = []
            # Convert source = {label_1: dict_1, ...} to results = [dict_1*, ...]
            # where dict_x* = dict_x + {label: label_x}
            for label, data_dict in source.items():
                result_dict = {'label': label}
                # Process all keys using get_value to handle H5MD format
                for key in data_dict.keys():
                    value = self.get_value(key, data_dict)
                    if value is not None:
                        result_dict[key] = value
                results.append(result_dict)
            return results

        return None


class H5MDArchiveWriter(MDParser):
    def __init__(self, **kwargs):
        self.h5_parser: H5MDH5Parser = H5MDH5Parser()
        self.simulation_parser = H5MDMetainfoParser()
        self.simulation_parser.max_nested_level = 10
        self.workflow_parser = H5MDMetainfoParser()
        super().__init__(**kwargs)

    def write_to_archive(self) -> None:
        # reload schema annotations
        reload(h5md)

        # create h5 parser
        self.h5_parser.filepath = self.mainfile

        trajectory_steps_data = Path(path='particles.all.position.step').get_data(
            self.h5_parser.data, default=[]
        )
        if trajectory_steps_data is None:
            trajectory_steps_data = []
        self.trajectory_steps = trajectory_steps_data
        self.h5_parser.trajectory_steps = self.trajectory_steps

        # TODO consider using a single parser for the whole archive
        # create metainfo parsers
        self.simulation_parser.annotation_key = h5md.HDF5_KEY
        simulation_data = Simulation()
        self.simulation_parser.data_object = simulation_data
        self.workflow_parser.annotation_key = h5md.HDF5_KEY
        workflow_data = molecular_dynamics.MolecularDynamics()
        self.workflow_parser.data_object = workflow_data

        # map from h5 source to metainfo target
        self.h5_parser.convert(self.simulation_parser)
        self.h5_parser.convert(self.workflow_parser)

        # assign simulation to archive data
        self.archive.data = self.simulation_parser.data_object
        self.archive.workflow2 = self.workflow_parser.data_object

        # close parsers
        self.h5_parser.close()
        self.simulation_parser.close()
        self.workflow_parser.close()

        # remove mapping annotations
        remove_mapping_annotations(self.archive.data.m_def)


class H5MDParser(MatchingParser):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.archive_writer = H5MDArchiveWriter()
        self.logger = get_logger(__name__)

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger | None = None,
        child_archives: dict[str, EntryArchive] | None = None,
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives or {})
