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
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import (
    MetainfoParser,
    Path,
    TextParser,
    XMLParser,
)
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import (
    SerialWorkflow,
    SimulationWorkflow,
    SinglePoint,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages.quantumespresso import common

from .common import libxc_shortcut, xc_functional_map
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


class XCFunctionalParser:
    @staticmethod
    def gen_string(data: dict[str, Any], separator='+') -> str:
        string = ''
        for key in sorted(data.keys()):
            val = data[key]
            weight = val.get('XC_functional_weight', 1.0)
            if string and weight > 0:
                string += separator
            if weight is not None:
                string += f'{weight:.3f}'
            string += val.get('XC_functional_name', '')
        return string

    @staticmethod
    def filter_data(data: dict[str, dict[str, Any]]) -> dict[str, Any]:
        out = dict()
        tol = 0.01
        for key, val in data.items():
            val_copy = val.copy()
            weight = val_copy.get('XC_functional_weight')
            if weight is None or abs(weight) < tol:
                continue
            else:
                if abs(weight - 1.0) < tol:
                    del val_copy['XC_functional_weight']
                val_copy.pop('exx_compute_weight', None)
            out[key] = val_copy
        return out


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

    def get_xc_functionals(self, source: str) -> list[dict[str, Any]]:
        numbers = source.split('(')[1].split(')')[0]
        nval = (4, 10)
        # handle different formatting
        if len(numbers) == nval[0]:
            # 4-digit format without spaces
            numbers_split = re.findall(r'(\d)', numbers)
        elif len(numbers) == nval[1]:
            # 5-digit format with/without spaces
            numbers_split = re.findall(r'[ \d]\d', numbers)
        else:
            # 6-digit with spaces
            numbers_split = numbers.split()

        if not numbers_split:
            self.logger.warning(
                'Unknown XC functional format', data=dict(value=numbers)
            )
            return []

        numbers_split = [int(n) for n in numbers_split]
        # numbers should have six digits
        numbers_split.extend([0] * (6 - len(numbers_split)))

        # map numbers to values
        xc_section_method = dict()
        xc_terms = dict()
        xc_terms_remove = dict()

        def get_data(source: list[dict[str, Any]]) -> dict[str, Any]:
            data = dict()
            exx_fraction = self.get_header('x_qe_exact_exchange_fraction', 0.0)
            for term in source:
                term_copy = term.copy()
                weight = term_copy.get('exx_compute_weight', 1.0)
                term_copy['XC_functional_weight'] = (
                    weight(exx_fraction) if not isinstance(weight, float) else weight
                )
                data.setdefault(term_copy.get('XC_functional_name', ''), term_copy)
            return data

        for i in range(6):
            xc_component = xc_functional_map[i]
            xc_number = numbers_split[i]
            if xc_number >= len(xc_component) or xc_component[xc_number] is None:
                continue
            xc_section_method.update(
                xc_component[xc_number].get('xc_section_method', {})
            )
            xc_terms.update(get_data(xc_component[xc_number].get('xc_terms', [])))
            xc_terms_remove.update(
                get_data(xc_component[xc_number].get('xc_terms_remove', []))
            )

        # remove terms
        for key, val in xc_terms_remove.items():
            weight = val.get('XC_functional_weight')
            xc_terms.setdefault(key, val)
            xc_terms[key]['XC_functional_weight'] *= -(weight or -1.0)

        # filter data
        xc_terms = XCFunctionalParser.filter_data(xc_terms)

        xc_functional_str = XCFunctionalParser.gen_string(xc_terms)
        if xc_functional_str in libxc_shortcut:
            # override for libXC compliance
            xc_terms = get_data(libxc_shortcut[xc_functional_str]['xc_terms'])
            xc_terms = XCFunctionalParser.filter_data(xc_terms)
            xc_functional_str = XCFunctionalParser.gen_string(xc_terms)
        # TODO make use of this
        xc_section_method['XC_functional'] = xc_functional_str

        return [xc_terms[key] for key in sorted(xc_terms.keys())]

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

    @property
    def program_name(self) -> str:
        return self.data_object.get('header', {}).get('program_name_version', [''])[0]

    @property
    def writers(self) -> Iterator[QuantumEspressoArchiveWriter]:
        if not self.data_object.get('program'):
            return

        for program in self.data.get('program', []):
            writer = load_writer(program[:30])
            if writer is None:
                self.logger.error('Parser not found for program.')
                continue
            writer.mainfile = self.filepath
            writer.mainfile_parser.data_object.mainfile = self.filepath
            # parse only the relevant program
            writer.mainfile_parser.data_object._file_handler = program.encode()
            yield writer


class MainfileXMLParser(XMLParser):
    _units_map = {'Hartree atomic units': dict(energy='hartee', length='bohr')}

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_datetime(self, date: str, time: str) -> datetime:
        return datetime.strptime(f'{date}{time}'.replace(' ', ''), '%d%b%Y%H:%M:%S')

    def apply_unit(self, value: np.ndarray | float, **kwargs) -> Any:
        unit = self._units_map.get(self.data.get('@Units')).get(kwargs.get('name'))
        if not unit or value is None:
            return value
        return value * ureg(unit)

    def get_forces(self, source: np.ndarray):
        return np.reshape(source, (np.size(source) // 3, 3))

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

    def parse_workflow(self) -> None:
        if self.archive.workflow2:
            return

        self.archive.workflow2 = SinglePoint()

        # multi run file
        multirun_workflow_archive = self.child_archives.get('workflow_multirun')
        if multirun_workflow_archive is not None:
            multirun_workflow_archive.workflow2 = SerialWorkflow(
                tasks=[TaskReference(task=self.archive.workflow2)]
            )
            for key, child_archive in self.child_archives.items():
                if key.startswith('workflow'):
                    continue
                multirun_workflow_archive.workflow2.tasks.append(
                    TaskReference(task=child_archive.workflow2)
                )

        # mainfiles in the same upload
        generic_workflow_archive = self.child_archives.get('workflow_generic')
        if (
            generic_workflow_archive is not None
            and self.archive.metadata.main_author is not None
        ):
            from nomad.app.v1.models import MetadataRequired  # noqa
            from nomad.search import search  # noqa

            parent_archive = multirun_workflow_archive or self.archive
            # add current archive workflow to generic workflow tasks
            generic_workflow_archive.workflow2 = SimulationWorkflow(
                tasks=[TaskReference(task=parent_archive.workflow2)]
            )

            upload_id = self.archive.metadata.upload_id
            metadata = search(
                owner='visible',
                user_id=self.archive.metadata.main_author.user_id,
                query={'upload_id': upload_id},
                required=MetadataRequired(
                    include=['entry_id', 'mainfile', 'parser_name']
                ),
            ).data
            parent_file = self.mainfile.split('raw/')[-1]
            parent_dir = os.path.dirname(parent_file)
            for result in metadata:
                parser_name = result.get('parser_name')
                # include only qe calculations
                if 'quantumespresso' not in parser_name:
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
                        == generic_workflow_archive.metadata.entry_id
                    ):
                        continue

                    # add workflow to generic workflow tasks
                    generic_workflow_archive.workflow2.tasks.append(
                        TaskReference(task=entry_archive.workflow2)
                    )
                    # add parent scf as input to task
                    if entry_archive.workflow2:
                        entry_archive.workflow2.inputs.append(
                            Link(section=parent_archive.workflow2)
                        )

    @property
    def mainfile_parser(self) -> MainfileTextParser | MainfileXMLParser:
        if self._mainfile_parser is None:
            ext = self.mainfile.rsplit('.', 1)[-1]
            self._mainfile_parser = dict(
                out=self._text_parser, xml=self._xml_parser
            ).get(ext)
            self._mainfile_parser.filepath = self.mainfile
            self.simulation_parser.annotation_key = dict(
                out=common.OUT_KEY, xml=common.XML_KEY
            ).get(ext)
        return self._mainfile_parser

    def write_to_archive(self) -> None:
        for n, writer in enumerate(self.mainfile_parser.writers):
            # write the first program to the main archive, the rest to child archives
            archive = (
                self.archive
                if n == 0
                else self.child_archives.get(
                    f'{n} {writer.mainfile_parser.program_name}'
                )
            )
            if archive is None:
                self.logger.error('Archive not found for program.')
                continue
            writer.parse_program(archive, n)

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

    archive_writer = QuantumEspressoArchiveWriter()
    _supported_exts = ['out', 'xml']

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
                        basename = os.path.basename(f).rsplit('.', 1)[1]
                        if basename not in basenames:
                            qe_files.append(f)
                if len(qe_files) > 1:
                    # generate workflow only if there is one scf file
                    qe_files = []
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
                children.extend(['workflow_multirun', *programs[1:]])

            self.creates_children = len(children) > 0

            return children or True

        return is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.level = len(child_archives)
        self.archive_writer.write(mainfile, archive, logger, child_archives)
