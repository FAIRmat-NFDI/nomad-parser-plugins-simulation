import os
from typing import TYPE_CHECKING, Iterable, Union

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from nomad.metainfo import (
        Context,
        Section,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

import nomad_simulations.schema_packages.properties.energies as energy_module
import nomad_simulations.schema_packages.properties.forces as force_module
import numpy as np

# from nomad.parsing.parser import MatchingParser
import structlog
from atomisticparsers.utils import MDParser
from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.file_parser import FileParser
from nomad.units import ureg
from nomad_simulations.schema_packages.general import (
    Program,
    Simulation,
)

# nomad-simulations
from nomad_simulations.schema_packages.general import Program as BaseProgram
from nomad_simulations.schema_packages.outputs import TotalEnergy, TotalForce

# # nomad-parser-gsd
# from nomad_parser_gsd.schema_packages.schema import (
#     Author,
#     EnergyEntry,
#     ForceEntry,
#     ModelSystem,
#     OutputsEntry,
#     ParamEntry,
#     Stress,
#     TrajectoryOutputs,
# )

logging = structlog.get_logger()
try:
    import gsd.fl as gsdfl
    import gsd.hoomd as gsdhoomd
    from gsd.hoomd import HOOMDTrajectory
    # import gsd.pygsd as gsdpy
except ImportError:
    logging.warn('Required module gsd.hoomd not found.')
    gsdhoomd = False
    # gsdpy = False

configuration = config.get_plugin_entry_point(
    'nomad_parser_gsd.parsers:parser_entry_point'
)

MOL = 6.022140857e23


class GSDFileParser(FileParser):
    def __init__(self):
        super().__init__(None)

    @property
    def filegsd(self):
        if self._file_handler is None:
            try:
                self._file_handler = gsdhoomd.open(name=self.mainfile, mode='r')
            except Exception:
                self.logger.error('Error reading gsd file.')

            if type(self._file_handler) is not HOOMDTrajectory:
                self.logger.error(
                    'Uknown GSD file object, only HOOMDTrajectory objects are supported.'
                )
        return self._file_handler

    def get_value(self, group, path: str, default=None):
        """
        Extracts group or dataset from group object based on path, and returns default if not defined.
        """
        section_segments = path.split('.')
        for section in section_segments:
            try:
                value = getattr(group, section)
                group = value
            except AttributeError:
                return

        return value if value is not None else default

    def parse(self, path: str = None, **kwargs):
        frame_path = '.'.join(path.split('.')[1:])
        frame = kwargs.get('frame', None)

        if not frame:
            value = None
        else:
            value = self.get_value(frame, frame_path)

        self._results[path] = value


class GSDParser(MDParser):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self._data_parser = GSDFileParser()
        self._maindir = None
        self._gsd_files = None
        self._basename = None
        self._program_dict = None
        self._n_frames = None
        self._time_unit = ureg.picosecond
        # ['N', 'position', 'orientation', 'types', 'typeid', 'mass', 'charge',
        # 'diameter', 'body', 'moment_inertia', 'velocity', 'angmom', 'image',
        # 'type_shapes']
        self._nomad_to_particles_group_map = {
            'n_atoms': 'N',
            'positions': 'position',
            'velocities': 'velocity',
            'forces': None,
            'labels': 'types',
            'label': None,
            'mass': 'mass',
            'charge': 'charge',
        }
        self._nomad_to_box_group_map = {
            'lattice_vectors': '_box',
            'periodic_boundary_conditions': None,
            'dimensionality': 'dimensions',
        }

    # Load GSD file as file layer object to access generating program name and version.
    def get_program_info(self):
        _file_layer = gsdfl.open(name=self.mainfile, mode='r')
        _application = _file_layer.application.split(' ')

        _program_info = dict()

        _program_info['name'] = _application[0]
        _program_info['version'] = _application[1]
        _program_info['schema'] = (
            _file_layer.schema
        )  # Name of the data schema. In example: 'hoomd'
        _program_info['schema_version'] = (
            _file_layer.schema_version
        )  # TODO: tuple(major, minor). If keep, format to major.minor?
        _program_info['gsd_version'] = (
            _file_layer.gsd_version
        )  # TODO: tuple(major, minor). Format?

        self._n_frames = _file_layer.nframes  # ? Not ideal to access in get_program, but only file layer object has attribute nframes

        _file_layer.close()

        return _program_info

    """
    Logged data encompasses values computed at simulation time that are too expensive 
    or cumbersome to re-compute in post processing. This specification does not define 
    specific chunk names or define logged data. Users may select any valid name for 
    logged data chunks as appropriate for their workflow.
    """

    # Additional data are stored in the log dictionary as numpy arrays:
    def get_logged_info(self):
        try:
            return gsdhoomd.read_log(name=self.mainfile, scalar_only=False)
        except FileNotFoundError:
            self.logger.warning(
                'No additional logged data found, no user-defined data will be stored.'
            )
            return dict()

    def get_system_info(self, frame_idx=None, frame=None):
        self._system_info = {'system': dict(), 'outputs': dict()}
        _path = f'{frame_idx}'
        _configuration_keys = ['step', 'dimensions', '_box']
        _interaction_keys = ['M', 'N', 'types', 'typeid', 'group', '_default_value']

        def get_particle_parameters(path: str = None, frame=None):
            n_particles = self._data_parser.get(f'{path}.N', frame=frame)
            if n_particles is None:
                return dict()
            else:
                return self._data_parser.get(path, frame=frame).__dict__

        def get_value(value, path=None, frame=None):
            if value is None:
                return value
            try:
                value = self._data_parser.get(
                    f'{path}.{value}' if path else value, frame=frame
                ).__dict__
                if value is None:
                    self.logger.warning(
                        f'No attributes found for key {value}. These attributes will '
                        'not be stored.'
                    )
                    return None
                else:
                    return value
            except KeyError:
                return value
            # TODO: handle data chunks that are not standard objects
            except AttributeError:
                pass

        _particle_data_dict = get_particle_parameters(
            path=f'{frame_idx}.particles', frame=frame
        )
        # print(_particle_data_dict.keys())
        if _particle_data_dict is None:
            self.logger.warning(
                f'No number of particles available in frame {frame_idx}. Other'
                ' particle attributes will not be stored for frame {frame_idx}.'
            )

        system_keys = {
            'time': ['system', 'outputs'],
            'n_atoms': 'system',
            'positions': 'system',
            'labels': 'system',
            'velocities': 'system',
            'mass': 'system',
            'charge': 'system',
            'lattice_vectors': 'system',
            'periodic_boundary_conditions': 'system',
            'dimensionality': 'system',
            'forces': 'outputs',
            'label': 'outputs',
        }

        # get quantities from particles chunk of GSD file
        for key, gsd_key in self._nomad_to_particles_group_map.items():
            section = system_keys[key]
            if isinstance(section, list):
                for sec in section:
                    self._system_info[sec][key] = (
                        _particle_data_dict[gsd_key] if gsd_key is not None else None
                    )
            else:
                self._system_info[system_keys[key]][key] = (
                    _particle_data_dict[gsd_key] if gsd_key is not None else None
                )
        # get steps and box quantities on from configurations chunk of GSD file:
        for key, gsd_key in self._nomad_to_box_group_map.items():
            section = system_keys[key]
            values_dict = get_value('configuration', path=_path, frame=frame)
            if isinstance(section, list):
                for sec in section:
                    self._system_info[sec][key] = (
                        values_dict[gsd_key] if gsd_key is not None else None
                    )
            else:
                self._system_info[system_keys[key]][key] = (
                    values_dict[gsd_key] if gsd_key is not None else None
                )

        # print(frame.__dict__.keys())
        # {'configuration': <gsd.hoomd.ConfigurationData object at 0x7f846452ad90>,
        # 'particles': <gsd.hoomd.ParticleData object at 0x7f846452ae10>,
        # 'bonds': <gsd.hoomd.BondData object at 0x7f846452add0>,
        # 'angles': <gsd.hoomd.BondData object at 0x7f846452ae90>,
        # 'dihedrals': <gsd.hoomd.BondData object at 0x7f846452af50>,
        # 'impropers': <gsd.hoomd.BondData object at 0x7f846452b010>,
        # 'constraints': <gsd.hoomd.ConstraintData object at 0x7f846452b0d0>,
        # 'pairs': <gsd.hoomd.BondData object at 0x7f846452b150>,
        # 'state': {},
        # 'log': {},
        # '_valid_state': ['hpmc/integrate/d', 'hpmc/integrate/a',
        # 'hpmc/sphere/radius', 'hpmc/sphere/orientable', 'hpmc/ellipsoid/a',
        # 'hpmc/ellipsoid/b', 'hpmc/ellipsoid/c', 'hpmc/convex_polyhedron/N',
        # 'hpmc/convex_polyhedron/vertices', 'hpmc/convex_spheropolyhedron/N',
        # 'hpmc/convex_spheropolyhedron/vertices',
        # 'hpmc/convex_spheropolyhedron/sweep_radius', 'hpmc/convex_polygon/N',
        # 'hpmc/convex_polygon/vertices', 'hpmc/convex_spheropolygon/N',
        # 'hpmc/convex_spheropolygon/vertices',
        # 'hpmc/convex_spheropolygon/sweep_radius',
        # 'hpmc/simple_polygon/N', 'hpmc/simple_polygon/vertices']}

        # for key in frame.__dict__.keys():
        #     values_dict = get_value(key, path=_path, frame=frame)
        #     if values_dict is not None:
        #         for value_key in values_dict.keys():
        #             print(value_key)

        _configuration_data_dict = get_value('configuration', path=_path, frame=frame)
        # self._system_info['system']['steps'] = _configuration_data_dict['step']
        # self._system_info['outputs']['steps'] = _configuration_data_dict['step']

        _bonds_dict = get_value('bonds', path=_path, frame=frame)
        _angles_dict = get_value('angles', path=_path, frame=frame)
        _dihedrals_dict = get_value('dihedrals', path=_path, frame=frame)
        _impropers_dict = get_value('impropers', path=_path, frame=frame)
        _pairs_dict = get_value('pairs', path=_path, frame=frame)

        return self._system_info

    def parse_system(self, simulation):
        system_info = self._system_info.get('system')
        if not system_info:
            self.logger.error('No particle information found in GSD file.')
            return

        # atoms_dict = system_info['step']
        # atoms_dict['is_representative'] = False

        # atom_labels = atoms_dict.get('labels')
        # if atom_labels is not None:
        #     try:
        #         # symbols2numbers(atom_labels)
        #         atoms_dict['labels'] = atom_labels
        #     except KeyError:  # TODO this check should be moved to the system normalizer in the new schema
        #         atoms_dict['labels'] = ['X'] * len(atom_labels)

        # topology = None
        # if i_step == 0:  # TODO extend to time-dependent bond lists and topologies
        #     atoms_dict['is_representative'] = True
        #     atoms_dict['bond_list'] = self._data_parser.get('connectivity.bonds')
        #     path_topology = 'connectivity.particles_group'
        #     topology = self._data_parser.get(path_topology)

        # # REMAP some of the data for the schema
        # atoms_dict['branch_label'] = (
        #     'Total System'  # ? Do we or should we have a default name for the entire system?
        # )
        # atoms_dict['time_step'] = atoms_dict.pop(
        #     'time'
        # ).magnitude  # TODO change in system_info
        # atomic_cell_keys = [
        #     'n_atoms',
        #     'lattice_vectors',
        #     'periodic_boundary_conditions',
        #     'positions',
        #     'velocities',
        #     'labels',
        # ]
        # atoms_dict['atomic_cell'] = {}
        # for key in atomic_cell_keys:
        #     atoms_dict['atomic_cell'][key] = atoms_dict.pop(key)

        # self.parse_trajectory_step(atoms_dict, simulation)

        # if i_step == 0 and topology:  # TODO extend to time-dependent topologies
        #     self.parse_system_hierarchy(
        #         simulation.model_system[-1], topology, path_topology
        #     )

    def write_to_archive(self) -> None:
        #######################################################################
        # Access simulation file(s).
        #######################################################################
        self._maindir = os.path.dirname(
            self.mainfile
        )  # ? GSD output single file or more?
        self._gsd_files = [
            _file for _file in os.listdir(self._maindir) if _file.endswith('.gsd')
        ]
        self._basename = os.path.basename(self.mainfile).rsplit('.', 1)[0]
        self._data_parser.mainfile = self.mainfile
        if self._data_parser.filegsd is None:
            self.logger.warning('GSD file missing in GSD Parser.')
            return

        #######################################################################
        # Start populating NOMAD-simulations schema
        #######################################################################
        simulation = Simulation()
        self._program_dict = self.get_program_info()
        simulation.program = Program(
            name=self._program_dict['name'],
            version=self._program_dict['version'],
        )
        for frame_idx, frame in enumerate(self._data_parser.filegsd):
            self.get_system_info(frame_idx=frame_idx, frame=frame)
            self.parse_system(simulation)
            print(self._system_info)
            # TODO: parse systems_info dict to nomad

            # ? Forces etc. are user-defined and read via get_logged_info?
            # TODO: Extract observables from logged data, parse to ModelOutput
            # system_keys = {
            #     'forces': 'calculation',
            # }
