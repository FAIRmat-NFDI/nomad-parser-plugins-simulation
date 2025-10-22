from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger
from nomad_simulation_parsers.parsers.utils.mdparserutils import MDParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.utils import get_logger
from nomad.units import ureg
from nomad.parsing.file_parser.mapping_parser import MappingParser, TextParser
from nomad_simulation_parsers.parsers.utils.mdanalysisparser import MDAnalysisParser
from typing import Any
from .log_parser import GromacsLogParser as GromacsLogTextParser
from .edr_parser import GromacsEDRParser as GromacsEDRFileParser
from nomad_simulations.schema_packages.general import Simulation
from importlib import reload
from nomad_simulation_parsers.schema_packages import gromacs
import os
import numpy as np
from ase.symbols import symbols2numbers

LOGGER = get_logger(__name__)

# TODO put this in utils
MOL = 6.022140857e23


class GromacsMetainfoParser(MetainfoParser):
    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GromacsLogParser(TextParser):
    _trajectory_steps_sampled: list[int] = []
    _trajectory_steps: list[int] = []
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


class GromacsEDRParser(MappingParser):
    # Fallback parser
    log_parser: GromacsLogParser = None
    edr_parser: GromacsEDRFileParser = None
    _trajectory_steps: list[int] = []
    _thermodynamic_steps: list[int] = []

    # TODO: temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self, **kwargs) -> dict[str | int, Any]:
        if self.data_object is not None:
            # Read data from log if cannot be parsed from edr
            for key in self.data_object.keys():
                self.data_object.parse(key)
            if not self.data_object.keys() and self.log_parser:
                data = {}
                steps = self.log_parser.input_parameters.get('step', [])
                for n, step in enumerate(steps):
                    energies = step.get('energies')
                    if energies is None:
                        continue
                    index = int(step.get('step_info', {}).get('Step', n))
                    for key in energies.keys():
                        data.setdefault(key, [None] * len(steps))
                        data[key][index] = step.energies.get(key)
                self.data_object._results = data
            return self.data_object._results
        return {}

    def from_dict(self, dct: dict[str, Any]):
        raise NotImplementedError

    def load_file(self) -> Any:
        if self.filepath:
            self.edr_parser.mainfile = self.filepath
        return self.edr_parser

    def get_outputs(self) -> list[dict[str, Any]]:
        outputs = []
        calculation_times = self.data.get('Time', [])
        energy_map = {
            'Total Energy': 'total',
            'Potential': 'potential',
            'Kinetic En.': 'kinetic',
            'pV': 'pressure_volume_work',
        }
        outputs = []
        for n, step in enumerate(self._thermodynamic_steps):
            data = dict(step=step, time=calculation_times[n] * ureg.picosecond)
            for key in self.data.keys():
                val = self.data.get(key)
                if val is None:
                    continue
                if key in energy_map:
                    val_n = val[n] * ureg.kilojoule / MOL
                    if key == 'Total Energy':
                        data.setdefault('energy', {})['value'] = val_n
                    else:
                        data.setdefault('energy', {}).setdefault(
                            'contributions', []
                        ).append(dict(value=val_n, label=energy_map[key]))
            if step in self._trajectory_steps:
                data['system_ref'] = f'/data/model_system/{self._trajectory_steps.index(step)}'
            outputs.append(data)
        return outputs


class GromacsMDAnalysisParser(MappingParser):
    aux_files: list[str] = []
    # parse pbc from log file
    log_parser: GromacsLogParser = None
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
            data = dict(
                forces=self.data_object.get_forces(n)
            )
            outputs.append(data)
        return outputs

    def get_configurations(self) -> list[dict[str, Any]]:
        configurations = []
        for step in self._trajectory_steps_sampled:
            n = self._trajectory_steps.index(step)
            data = dict(
                labels=self.get_atom_labels(n),
                positions=self.data_object.get_positions(n),
                velocities=self.data_object.get_velocities(n),
                lattice_vectors=self.data_object.get_lattice_vectors(n),
            )
            if self.log_parser is not None:
                pbc = self.log_parser.input_parameters.get('pbc', 'xyz')
                data['pbc'] = ['x' in pbc, 'y' in pbc, 'z' in pbc]
            configurations.append(data)
        return configurations


class GromacsArchiveWriter(MDParser):
    def __init__(self, **kwargs):
        self._simulation_parser = GromacsMetainfoParser()
        self._log_parser = GromacsLogParser(text_parser=GromacsLogTextParser())
        self._edr_parser = GromacsEDRParser(edr_parser=GromacsEDRFileParser())
        self._mdanalysis_parser = GromacsMDAnalysisParser(mdanalysis_parser=MDAnalysisParser())
        super().__init__(**kwargs)

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
        self._mdanalysis_parser._trajectory_steps_sampled = self.trajectory_steps
        self._mdanalysis_parser._trajectory_steps = traj_steps
        self._mdanalysis_parser._thermodynamic_steps = self.thermodynamics_steps
        self._edr_parser._thermodynamic_steps = self.thermodynamics_steps
        self._edr_parser._trajectory_steps = traj_steps

        # parse main log file
        self._simulation_parser.annotation_key = 'log'
        self._log_parser.convert(self._simulation_parser)


        # parse edr file
        self._simulation_parser.annotation_key = 'edr'
        self._edr_parser.log_parser = self._log_parser
        self._edr_parser.convert(self._simulation_parser)

        # parse mdanalysis trajectory files
        self._simulation_parser.annotation_key = 'tpr'
        self._mdanalysis_parser.log_parser = self._log_parser
        self._mdanalysis_parser.convert(self._simulation_parser)

        self._simulation_parser.close()
        self._log_parser.close()
        self._mdanalysis_parser.close()
        self._edr_parser.close()

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
