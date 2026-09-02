from __future__ import annotations

import os
import re
from collections.abc import Iterable, Iterator
from datetime import datetime
from types import ModuleType
from typing import Any

import numpy as np
from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Link, TaskReference
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import (
    MetainfoParser,
    Path,
    TextParser,
    XMLParser,
)
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import (
    SerialWorkflow,
    SimulationWorkflow,
    SinglePoint,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import (
    link_outputs_to_model_systems,
    search_files,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import common

from .file_parser import QuantumEspressoFileParser

LOGGER = get_logger(__name__)
PROGRAM_NAME_RE = re.compile(
    r'(?:Program +(\w+) +v\.([\d\.]+))|(?:<creator NAME=\"(\w+?)\" VERSION=\"([\d\.]+))'
)


# TODO temporary fix for structlog unable to propagate logger
class QuantumEspressoMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


def get_program_name_version(header: str) -> tuple[str, tuple[int]]:
    match = PROGRAM_NAME_RE.search(header)
    if not match:
        return ('', ())
    name = (match.group(1) or match.group(3)).lower()
    version = tuple(
        [int(v) for v in (match.group(2) or match.group(4)).split('.') if v.isdecimal()]
    )
    return name, version


def load_writer(header: str) -> QuantumEspressoArchiveWriter:
    from .epw.parser import EPWArchiveWriter  # noqa
    from .gipaw.parser import GIPAWArchiveWriter  # noqa
    from .phonon.parser import PhononArchiveWriter  # noqa
    from .pwscf.parser import PWSCFArchiveWriter  # noqa
    from .xspectra.parser import XSpectraArchiveWriter  # noqa

    _writers = {
        'pwscf': PWSCFArchiveWriter(),
        'epw': EPWArchiveWriter(),
        'phonon': PhononArchiveWriter(),
        'xspectra': XSpectraArchiveWriter(),
        'gipaw': GIPAWArchiveWriter(),
    }
    name = get_program_name_version(header)
    return _writers.get(name[0] if name[0] else header)


class MainfileTextParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_version(self, name_version: list[str]) -> str:
        return ' '.join(name_version[1:]).lstrip('v.')

    def get_datetime(self, date_time: str) -> datetime:
        return datetime.strptime(date_time.replace(' ', ''), '%d%b%Y%H:%M:%S')

    def get_header(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, self.data.get('header', {}).get(key, default))

    def get_n_spin_channels(self):
        magnetic = self.get_header('starting_magnetization')
        if magnetic is None:
            calculation = self.data.get('self_consistent', {})
            magnetic = calculation.get(
                'magnetization_total', calculation.get('spin_pol')
            )
        return 1 if magnetic is None else 2

    def get_energy_contributions(self, source: dict[str, Any]):
        return [
            dict(value=val.magnitude, name=key.split('energy_total_', 1)[-1])
            for key, val in source.items()
            if key != 'energy_total' and val is not None
        ]

    # Quantum ESPRESSO prints the XC either as a single functional name ('PBE')
    # or as the four DFT slot codes ('SLA PW PBX PBC'), optionally followed by
    # the numeric slot codes in parentheses. Map the slot combination to a
    # standard functional name; the schema expands it into LibXC components and
    # derives `jacobs_ladder`.
    _slot_combo_names = {
        'SLA PZ': 'PZ81',
        'SLA PW': 'LDA',
        'SLA VWN': 'VWN',
        'SLA PW PBX PBC': 'PBE',
        'SLA PW PSX PSC': 'PBEsol',
        'SLA PW HHNX PBC': 'RPBE',
        'SLA PW REVX PBC': 'revPBE',
        'SLA PW WCX PBC': 'WC',
        'SLA PW GGX GGC': 'PW91',
        'SLA PW B88 P86': 'BP86',
        'SLA LYP B88 BLYP': 'BLYP',
    }

    def get_functional_key(self, source: str) -> str | None:
        if not source:
            return None
        # Drop the trailing numeric slot codes '( 1 4 3 4 0 0)' when present.
        tokens = source.split('(', 1)[0].split()
        if not tokens:
            return None
        if len(tokens) == 1:
            return tokens[0]
        combo = ' '.join(tokens)
        name = self._slot_combo_names.get(combo)
        if name is None:
            # preserve the raw slot combination so the reported XC is not dropped
            self.logger.debug('unmapped QE XC slot combination', data=dict(value=combo))
            return combo
        return name

    def get_value(self, source: dict[str, Any], key: str = '', units: str = 'units'):
        key_split = key.rsplit('.', 1)
        parent = Path(path=key_split[0]).get_data(source)
        header = self.data.get('header', {})
        if parent is None:
            source = header
            parent = self.get_value(header, key_split[0], '')
        value = parent if len(key_split) == 1 else parent.get(key_split[1])

        if value is None or not units:
            return list(value) if isinstance(value, np.ndarray) else value

        units = (source if len(key_split) == 1 else parent).get(units, units).lower()
        alat = source.get('alat', header.get('alat', 1.0))
        value = np.array(value, dtype=float)

        if units in ['alat', 'a_0']:
            value *= alat
        elif units in ['bohr', 'angstrom']:
            units_mapping = dict(bohr=ureg.bohr, angstrom=ureg.angstrom)
            value = value * units_mapping.get(units)
        elif units == '2 pi/alat':
            value *= 2 * np.pi / alat
        elif units == 'crystal':
            cell = self.get_value(source, 'simulation_cell', '')
            if cell is not None:
                value = np.dot(
                    value.magnitude if hasattr(value, 'magnitude') else value,
                    cell.magnitude if hasattr(cell, 'magnitude') else cell,
                ) * getattr(cell, 'units', 1.0)
        return value

    def get_periodic_boundary_conditions(
        self, source: dict[str, Any]
    ) -> list[bool] | None:
        cell = self.get_value(source, 'simulation_cell', '')
        if cell is None:
            return None
        return [True, True, True]

    @property
    def program_name(self) -> str:
        return self.data_object.get('header', {}).get('program_name_version', [''])[0]

    @property
    def writers(self) -> Iterator[QuantumEspressoArchiveWriter]:
        if not self.data_object.get('program'):
            return

        start = 0
        for program in self.data.get('program', []):
            writer = load_writer(program.header[:30])
            if writer is None:
                self.logger.error('Parser not found for program.')
                continue
            writer.mainfile = self.filepath
            pointers = program._file_handler[0]
            if isinstance(writer.mainfile_parser, TextParser):
                writer.mainfile_parser.data_object.mainfile = self.filepath
                # parse only the relevant program
                writer.mainfile_parser.data_object._file_handler = [
                    (pointers[0] + start, pointers[1] + start)
                ]
            start += pointers[1]
            yield writer


class MainfileXMLParser(XMLParser):
    _units_map = {'Hartree atomic units': dict(energy='hartree', length='bohr')}

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_datetime(self, date: str, time: str) -> datetime:
        return datetime.strptime(f'{date}{time}'.replace(' ', ''), '%d%b%Y%H:%M:%S')

    def apply_unit(self, value: np.ndarray | float, **kwargs) -> Any:
        unit = self._units_map.get(self.data.get('@Units'), {}).get(kwargs.get('name'))
        if not unit or value is None:
            return value
        return value * ureg(unit)

    def get_forces(self, source: np.ndarray):
        return np.reshape(source, (np.size(source) // 3, 3))

    def get_periodic_boundary_conditions(self, cell: Any = None) -> list[bool] | None:
        if cell is None:
            return None
        return [True, True, True]

    def get_energy_contributions(self, source: dict[str, Any]):
        return [
            dict(value=val, name=key) for key, val in source.items() if key != 'etot'
        ]

    @property
    def program_name(self) -> str:
        keys = list(self.data.keys())
        source = self.data if 'general_info' in keys else self.data[keys[0]]
        return source.get('general_info', {}).get('creator', {}).get('@NAME', '')

    @property
    def writers(self) -> Iterator[QuantumEspressoArchiveWriter]:
        for dct in self.data.values():
            writer = load_writer(self.program_name.lower())
            if writer is None:
                self.logger.error('Parser not found for program.')
                continue
            writer.mainfile = self.filepath
            writer.mainfile_parser._data = dct
            yield writer


class QuantumEspressoArchiveWriter(ArchiveWriter):
    """
    Wrapper for the program-specific archive writer.
    """

    schema: ModuleType = common
    simulation_parser = QuantumEspressoMetainfoParser()
    _text_parser = MainfileTextParser(text_parser=QuantumEspressoFileParser())
    _xml_parser = MainfileXMLParser()
    _mainfile_parser = None

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.data_object = Simulation(
            program=Program(name='Quantum Espresso')
        )
        # convert
        self.mainfile_parser.convert(self.simulation_parser)
        # set the parsed data to archive
        archive.data = self.simulation_parser.data_object

        if not archive.workflow2:
            # set workflow to single point by default
            archive.workflow2 = SinglePoint()

    def _link_modules(self) -> None:
        """
        Link archives created from multiple modules in one mainfile.
        """

        workflow_archive = self.child_archives.get('workflow_modules')
        if not workflow_archive:
            return

        workflow_archive.workflow2 = SerialWorkflow(
            tasks=[TaskReference(task=self.archive.workflow2)]
        )
        for key, child_archive in self.child_archives.items():
            if key.startswith('workflow'):
                continue
            workflow_archive.workflow2.tasks.append(
                TaskReference(task=child_archive.workflow2)
            )

    def _link_files(self) -> None:
        """
        Link archives created from separate mainfiles in the same upload.
        """
        workflow_archive = self.child_archives.get('workflow_generic')
        if workflow_archive is None or self.archive.metadata.main_author is None:
            return

        from nomad.app.v1.models import MetadataRequired  # noqa
        from nomad.search import search  # noqa

        parent_archive = self.child_archives.get('workflow_modules') or self.archive
        # add current archive workflow to generic workflow tasks
        workflow_archive.workflow2 = SimulationWorkflow(
            tasks=[TaskReference(task=parent_archive.workflow2)]
        )

        upload_id = self.archive.metadata.upload_id
        metadata = search(
            owner='visible',
            user_id=self.archive.metadata.main_author.user_id,
            query={'upload_id': upload_id},
            required=MetadataRequired(include=['entry_id', 'mainfile', 'parser_name']),
        ).data
        parent_file = self.mainfile.split('raw/')[-1]
        parent_dir = os.path.dirname(parent_file)
        for result in metadata:
            # include only qe calculations
            if 'quantumespresso' not in result.get('parser_name'):
                continue

            entry_id = result.get('entry_id')
            if not entry_id or entry_id == self.archive.metadata.entry_id:
                # skip the current entry
                continue

            # link only entries in the same directory or sub-directories
            mainfile = result.get('mainfile')
            if mainfile and mainfile.startswith(parent_dir):
                entry_archive: EntryArchive = self.archive.m_context.load_archive(
                    entry_id, upload_id, None
                )
                if not entry_archive.workflow2:
                    continue
                if (
                    entry_archive.metadata.entry_id
                    == workflow_archive.metadata.entry_id
                ):
                    continue

                if entry_archive.data and entry_archive.data.outputs:
                    # add model_system refs
                    if self.archive.data.model_system:
                        entry_archive.data.outputs[
                            0
                        ].model_system_ref = self.archive.data.model_system[0]

                # add workflow to generic workflow tasks
                workflow_archive.workflow2.tasks.append(
                    TaskReference(task=entry_archive.workflow2)
                )
                # add parent scf as input to task
                if entry_archive.workflow2:
                    entry_archive.workflow2.inputs.append(
                        Link(section=parent_archive.workflow2)
                    )

    def parse_workflow(self) -> None:
        self._link_modules()
        self._link_files()

    @property
    def mainfile_parser(self) -> MainfileTextParser | MainfileXMLParser:
        if self._mainfile_parser is None:
            ext = self.mainfile.rsplit('.', 1)[-1].lower()
            self._mainfile_parser = dict(
                out=self._text_parser, log=self._text_parser, xml=self._xml_parser
            ).get(ext)
            if self._mainfile_parser is None:
                self.logger.error('Parser not found for mainfile extension.')
                return None
            self._mainfile_parser.filepath = self.mainfile
            self.simulation_parser.annotation_key = dict(
                out=common.OUT_KEY, log=common.OUT_KEY, xml=common.XML_KEY
            ).get(ext)
        return self._mainfile_parser

    def write_to_archive(self) -> None:
        for n, writer in enumerate(self.mainfile_parser.writers):
            # write the first program to the main archive, the rest to child archives
            archive = (
                self.archive
                if n == 0
                else self.child_archives.get(
                    f'{n} {writer.mainfile_parser.program_name.lower()}'
                )
            )
            if archive is None:
                self.logger.error('Archive not found for program.')
                continue
            writer.parse_program(archive, n)
            writer.mainfile_parser.close()

        self.mainfile_parser.close()
        self.simulation_parser.close()

        self.parse_workflow()


def sort_qe_files(filenames: list[str]) -> list[tuple[str, datetime]]:
    """
    Sort QE mainfiles based on execution time.
    """
    sorted_files = []
    re_pattern = re.compile(
        r'starts on *(\w+) *at *([\d ]+\:[\d ]+\:[\d ]+)|'
        r'DATE\="(\w+)"\s+TIME\="([\d ]+\:[\d ]+\:[\d ]+)"'
    )
    for name in filenames:
        with open(name) as f:
            head = f.read(config.process.parser_matching_size)
            match = re_pattern.search(head)
            if not match:
                continue
            sorted_files.append(
                (
                    name,
                    datetime.strptime(
                        ''.join([g for g in match.groups() if g]).replace(' ', ''),
                        '%d%b%Y%H:%M:%S',
                    ),
                )
            )

    sorted_files.sort(key=lambda x: x[1])
    return sorted_files


def get_program_types(filename: str, multiple=False) -> list[str]:
    """
    Determine type of program(s) in file. If multiple will read all programs.
    """
    programs = []
    with open(filename) as f:
        for line in f:
            name_version = get_program_name_version(line)
            if name_version[0]:
                programs.append(name_version[0])
                if not multiple:
                    break
    return programs


class QuantumEspressoParser(MatchingParser):
    """
    Common parser for Quantum Espresso mainfiles including
    PWSCF, Phonon, EPW and XSpectra.
    """

    _supported_exts = ['out', 'log', 'xml']

    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable[str]:
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
        if is_mainfile:
            children = []
            programs = get_program_types(filename, multiple=True)
            if not programs:
                return True
            if 'pwscf' in programs[0].lower():
                # search all qe mainfiles in the directory and sub directories
                qe_files = []
                basenames = []
                for ext in self._supported_exts:
                    for f in search_files(
                        f'*.{ext}', os.path.dirname(filename), include_all=True
                    ):
                        basename = os.path.basename(f).rsplit('.', 1)[0]
                        if basename not in basenames:
                            basenames.append(basename)
                            qe_files.append(f)
                if len(qe_files) > 1:
                    # generate workflow only if there is one scf file
                    other_programs = []
                    for f in qe_files:
                        other_programs.extend([p.lower() for p in get_program_types(f)])
                    if other_programs.count('pwscf') > 1:
                        LOGGER.warning(
                            """Found multiple PWSCF files. Not generating workflow"""
                        )
                    else:
                        children.append('workflow_generic')

            if len(programs) > 1:
                # create separate entries for each program instance
                children.extend(
                    [
                        'workflow_modules',
                        *[f'{n + 1} {name}' for n, name in enumerate(programs[1:])],
                    ]
                )

            self.creates_children = len(children) > 0

            return children or True

        return is_mainfile

    def parse(  # noqa: PLR0912
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.level = len(child_archives)
        archive_writer = QuantumEspressoArchiveWriter()
        archive_writer.write(mainfile, archive, logger, child_archives)

        # TODO add this in the archive writer
        if archive.data and archive.data.outputs:
            link_outputs_to_model_systems(archive.data)
