import os
from importlib import reload
from typing import Any

import numpy as np
from ase.symbols import symbols2numbers
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser.mapping_parser import (
    MappingParser,
    MetainfoParser,
    TextParser,
)
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser
from nomad_simulation_parsers.schema_packages import gromacs

from .edr_parser import GromacsEDRParser as GromacsEDRFileParser
from .log_parser import GromacsLogParser as GromacsLogTextParser
from .mdanalysis_parser import GromacsMDAnalysisParser as MDAnalysisParser
from .mdp_parser import GromacsMdpParser as GromacsMDPTextParser
from .xvg_parser import GromacsXvgParser as GromacsXVGTextParser

LOGGER = get_logger(__name__)
ENERGY_UNIT = ureg.kilojoule / ureg.avogadro_number


class GromacsMetainfoParser(MetainfoParser):
    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GromacsThermodynamicsParser(MappingParser):
    _trajectory_steps: list[int] = []
    _thermodynamic_steps: list[int] = []
    _base_calc_unit_map = {
        'Temperature': ureg.kelvin,
        'Volume': ureg.nm**3,
        'Density': ureg.kilogram / ureg.m**3,
        'Pressure (bar)': ureg.bar,
        'Pressure': ureg.bar,
        'Enthalpy': ENERGY_UNIT,
    }
    _energy_label_map = {
        'Potential': 'potential',
        'Kinetic En.': 'kinetic',
        'Total Energy': 'total',
        'pV': 'pressure_volume_work',
    }

    def get_outputs(self, source: dict[str, Any] = {}) -> list[dict[str, Any]]:
        outputs = []
        times = source.get('Time', [])
        outputs = []
        for n, step in enumerate(self._thermodynamic_steps):
            data = dict(
                step=step,
                time=times[n] * ureg.picosecond if times[n] is not None else None,
            )
            for key in source.keys():
                val = source.get(key)
                if val is None or val[n] is None:
                    continue
                if key in self._energy_label_map:
                    val_n = val[n] * ENERGY_UNIT
                    if key == 'Total Energy':
                        data.setdefault('energy', {})['value'] = val_n
                    else:
                        data.setdefault('energy', {}).setdefault(
                            'contributions', []
                        ).append(dict(value=val_n, label=self._energy_label_map[key]))
                elif key in self._base_calc_unit_map:
                    data[key] = val[n] * self._base_calc_unit_map[key]
            if step in self._trajectory_steps:
                data['system_ref'] = (
                    f'/data/model_system/{self._trajectory_steps.index(step)}'
                )
            outputs.append(data)
        return outputs


class GromacsLogParser(TextParser, GromacsThermodynamicsParser):
    _trajectory_steps_sampled: list[int] = []
    input_parameters: dict[str, Any] = {}

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self):
        data_object = super().load_file()

        def normalize(input_dict: dict[str, Any]) -> None:
            for key, val in input_dict.items():
                if isinstance(val, dict):
                    normalize(val)
                elif isinstance(val, str):
                    input_dict[key.replace('_', '-')] = val.lower()
                elif isinstance(val, float):
                    if abs(val) == np.inf:
                        input_dict[key] = 'inf' if val > 0 else '-inf'

        self.input_parameters = data_object.get('input_parameters', {})
        normalize(self.input_parameters)

        return data_object

    def get_version(self, source: str) -> str:
        return (source or 'unknown').lstrip('VERSION')

    def get_configurations(self):
        pbc = [k in self.input_parameters.get('pbc', 'xyz') for k in 'xyz']
        return [dict(pbc=pbc) for _ in self._trajectory_steps_sampled]

    def get_outputs(self) -> list[dict[str, Any]]:
        data = {}
        steps = self.data.get('step', [])
        n_steps = len(self._thermodynamic_steps)
        for n, step in enumerate(steps):
            energies = step.get('energies')
            info = step.get('step_info')
            if energies is None or info is None:
                continue
            step_n = int(info.get('Step'))
            if step_n not in self._thermodynamic_steps:
                continue
            index = self._thermodynamic_steps.index(step_n)
            for key in energies.keys():
                data.setdefault(key, [None] * n_steps)
                data[key][index] = energies.get(key)
            data.setdefault('Time', [None] * n_steps)
            data['Time'][index] = info.get('Time', 1)
        outputs = super().get_outputs(data)
        return outputs


class GromacsMDPParser(TextParser):
    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self):
        data_object = super().load_file()

        def check_input_parameters_dict_recursive(self, input_dict: dict, key: str):
            if key in input_dict:
                return True
            for _, v in input_dict.items():
                if isinstance(v, dict):
                    if check_input_parameters_dict_recursive(v, key):
                        return True
            return False

        input_parameters = data_object.get('input_parameters', {})
        # parameters that are unique to the mdp file
        input_parameters['mdp_unique_params'] = {}
        for key, param in input_parameters.items():
            new_key = key.replace('_', '-')
            if not check_input_parameters_dict_recursive(input_parameters, new_key):
                input_parameters['mdp_unique_params'][new_key] = (
                    param.lower() if isinstance(param, str) else param
                )

        self.input_parameters = input_parameters
        return data_object


class GromacsEDRParser(GromacsThermodynamicsParser):
    # Fallback parser
    edr_parser: GromacsEDRFileParser = None

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs) -> dict[str | int, Any]:
        if self.data_object is not None:
            for key in self.data_object.keys():
                self.data_object.parse(key)
        return self.data_object._results

    def from_dict(self, dct: dict[str, Any]):
        raise NotImplementedError

    def load_file(self) -> Any:
        if self.filepath:
            self.edr_parser.mainfile = self.filepath
        return self.edr_parser


class GromacsXVGParser(TextParser):
    input_parameters = {}

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_results(self) -> dict[str, Any]:
        title = self.data.get('title', '')
        # TODO incorporate x and y axis labels into the checks
        if not (r'dH/d\xl\f{}' in title and r'\xD\f{}H' in title):
            return {}
        results = {}

        free_energy = results.setdefault('free_energy_calculations', {})
        columns = self.data.get('column_vals')
        free_energy['n_frames'] = len(columns)
        free_energy['value_unit'] = str(ENERGY_UNIT.units)
        # TODO get n_states from input_parameters
        free_energy['n_states'] = self.input_parameters.get('n_states')
        xaxis = self.data.get('xaxis', '').lower()
        # The expected columns of the xvg file are:
        # Total Energy
        # dH/dlambda current lambda
        # Delta H between each lambda and current lambda (n_lambda columns)
        # PV Energy
        if 'time' in xaxis and columns[:, 3:-1].shape[1] == free_energy['n_states']:
            free_energy['times'] = columns[:, 0] * ureg.ps
            columns = columns[:, 1:] * ENERGY_UNIT.magnitude

        return results


class GromacsMDAnalysisParser(MappingParser):
    aux_files: list[str] = []
    mdanalysis_parser: MDAnalysisParser = None
    _trajectory_steps_sampled: list[int] = []
    _trajectory_steps: list[int] = []
    _thermodynamic_steps: list[int] = []

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs) -> dict[str | int, Any]:
        if self.data_object is not None:
            self.data_object.parse()
            return self.data_object._results
        return {}

    def from_dict(self, dct: dict[str, Any]):
        raise NotImplementedError

    def load_file(self) -> MDAnalysisParser:
        if self.filepath:
            self.mdanalysis_parser.mainfile = self.filepath
            if self.aux_files:
                self.mdanalysis_parser.auxilliary_files = self.aux_files
        return self.mdanalysis_parser

    def get_atom_labels(self, index: int = 0) -> list[str]:
        labels = self.data_object.get_atom_labels(index)
        try:
            symbols2numbers(labels)
        except Exception:
            labels = ['X'] * len(labels)
        return labels

    def get_outputs(self) -> list[dict[str, Any]]:
        outputs = []
        for step in self._thermodynamic_steps:
            if step not in self._trajectory_steps:
                continue
            n = self._trajectory_steps.index(step)
            data = dict(forces=self.data_object.get_forces(n))
            outputs.append(data)
        return outputs

    def get_configurations(self) -> list[dict[str, Any]]:
        configurations = []
        for n, _ in enumerate(self._trajectory_steps_sampled):
            configurations.append(
                dict(
                    labels=self.get_atom_labels(n),
                    positions=self.data_object.get_positions(n),
                    velocities=self.data_object.get_velocities(n),
                    lattice_vectors=self.data_object.get_lattice_vectors(n),
                )
            )
        return configurations


class GromacsArchiveWriter(MDParser):
    def __init__(self, **kwargs):
        self._simulation_parser = GromacsMetainfoParser()
        self._log_parser = GromacsLogParser(text_parser=GromacsLogTextParser())
        self._mdp_parser = GromacsMDPParser(text_parser=GromacsMDPTextParser())
        self._edr_parser = GromacsEDRParser(edr_parser=GromacsEDRFileParser())
        self._mdanalysis_parser = GromacsMDAnalysisParser(
            mdanalysis_parser=MDAnalysisParser()
        )
        self._xvg_parser = GromacsXVGParser(text_parser=GromacsXVGTextParser())
        self.mdp_ext = 'mdp'
        self.mdp_std_filename = 'mdout'
        super().__init__(**kwargs)

    def get_mdp_file(self):
        """
        Tries to find the mdp input parameters (ext = mdp) that match the mainfile
        calculation.
        Priority is as follows:
            1. output mdp file containing both the matching mainfile name and the
            standard gromacs name `mdout`
            2. file containing the standard gromacs name `mdout`
            3. input mdp file matching the mainfile name (as usual)
            4. any `.mdp` file within the directory (as usual)
        """
        files = [d for d in self._gromacs_files if d.endswith(self.mdp_ext)]

        if len(files) == 0:
            return ''

        if len(files) == 1:
            return os.path.join(self._maindir, files[0])

        for f in files:
            filename = f.rsplit('.', 1)[0]
            if self._basename in filename and self.mdp_std_filename in filename:
                return os.path.join(self._maindir, f)

        for f in files:
            filename = f.rsplit('.', 1)[0]
            if self.mdp_std_filename in filename:
                return os.path.join(self._maindir, f)

        return self.get_gromacs_file(self.mdp_ext)

    def get_gromacs_file(self, ext: str) -> str:
        files = [d for d in self._gromacs_files if d.endswith(ext)]

        if len(files) == 0:
            return ''

        if len(files) == 1:
            return os.path.join(self._maindir, files[0])

        basenames = [f.rsplit('.', 1)[0] for f in files]
        # we assume that the file has the same basename as the log file e.g.
        # out.log would correspond to out.tpr and out.trr and out.edr
        for n, basename in enumerate(basenames):
            if basename == self._basename:
                return os.path.join(self._maindir, files[n])

        for n, basename in enumerate(basenames):
            if basename.startswith(self._basename):
                return os.path.join(self._maindir, files[n])

        # if the files are all named differently, we guess that the one that does not
        # share the same basename would be file we are interested in
        # e.g. out.log someout.log out.tpr out.trr another.tpr file.trr
        # we guess that the out.* files belong together and the rest that does not share
        # a basename would be grouped together
        counts = []
        all_basenames = [f.rsplit('.', 1)[0] for f in self._gromacs_files]
        for n, basename in enumerate(basenames):
            count = 0
            for ref_basename in all_basenames:
                if basename == ref_basename:
                    count += 1
            if count == 1:
                return os.path.join(self._maindir, files[n])
            counts.append(count)

        return os.path.join(self._maindir, files[counts.index(min(counts))])

    def write_to_archive(self):
        # intitialize variables
        self._maindir = os.path.dirname(self.mainfile)
        self._gromacs_files = os.listdir(self._maindir)
        self._basename = os.path.basename(self.mainfile).rsplit('.', 1)[0]
        # reload the schema annotations
        reload(gromacs)

        self.archive.data = Simulation(program=Program(name='GROMACS'))
        self._simulation_parser.data_object = self.archive.data

        # set up source parsers
        self._log_parser.filepath = self.mainfile
        self._edr_parser.filepath = self.get_gromacs_file('edr')
        self._mdanalysis_parser.filepath = self.get_gromacs_file('tpr')
        # TODO include input parameters read from mdp parser
        self._mdp_parser.filepath = self.get_mdp_file()
        self._xvg_parser.filepath = self.get_gromacs_file('xvg')
        # determine auxiliary file order: trr, xtc
        for ext in ['trr', 'xtc']:
            aux_file = self.get_gromacs_file(ext)
            if not aux_file:
                continue
            self._mdanalysis_parser.aux_files = [aux_file]
            # check if positions are parsed
            positions = self._mdanalysis_parser.data_object.get_positions(0)
            if positions is not None:
                break

        # determine sampled trajectory steps
        n_frames = self._mdanalysis_parser.data_object.get('n_frames', 0)
        traj_sampling_rate = self._log_parser.input_parameters.get('nstxout', 1)
        self.n_atoms = [
            self._mdanalysis_parser.data_object.get_n_atoms(n) for n in range(n_frames)
        ]
        traj_steps = [n * traj_sampling_rate for n in range(n_frames)]
        self.trajectory_steps = traj_steps

        # determine sampled thermodynamic steps
        calculation_times = self._edr_parser.data.get('Time', [])
        time_step = self._log_parser.input_parameters.get('dt')
        if time_step is None and len(calculation_times) > 1:
            time_step = calculation_times[1] - calculation_times[0]
        self.thermodynamics_steps = [
            int(time / time_step if time_step else 1) for time in calculation_times
        ]

        # pass sampled steps to parsers
        self._log_parser._trajectory_steps_sampled = self.trajectory_steps
        self._log_parser._trajectory_steps = traj_steps
        self._log_parser._thermodynamic_steps = self.thermodynamics_steps
        self._mdanalysis_parser._trajectory_steps_sampled = self.trajectory_steps
        self._mdanalysis_parser._trajectory_steps = traj_steps
        self._mdanalysis_parser._thermodynamic_steps = self.thermodynamics_steps
        self._edr_parser._thermodynamic_steps = self.thermodynamics_steps
        self._edr_parser._trajectory_steps = traj_steps

        # parse main log file
        self._simulation_parser.annotation_key = gromacs.LOG_KEY
        self._log_parser.convert(self._simulation_parser)

        # parse edr file
        self._simulation_parser.annotation_key = gromacs.EDR_KEY
        self._edr_parser.convert(self._simulation_parser)

        # parse mdanalysis trajectory files
        self._simulation_parser.annotation_key = gromacs.TPR_KEY
        self._mdanalysis_parser.convert(self._simulation_parser)

        self._simulation_parser.close()
        self._log_parser.close()
        self._mdanalysis_parser.close()
        self._edr_parser.close()
        self._mdp_parser.close()
        self._xvg_parser.close()


class GromacsParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = GromacsArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = None,
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
