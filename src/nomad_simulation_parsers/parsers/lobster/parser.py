import datetime
import os
from typing import Any

import ase
import ase.io
import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.file_parser import FileParser
from nomad_file_parser.mapping_parser import (
    MappingParser,
    MetainfoParser,
    TextParser,
)
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import lobster

from .file_parser import (
    CHARGEParser,
    COXPCARParser,
    DOSCARParser,
    ICOXPLISTParser,
    OutParser,
)

LOGGER = get_logger(__name__)


class LobsterMainfileParser(TextParser):
    text_parser = OutParser()

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_unix_time(self, datetime_str: str):
        return (
            datetime.datetime.strptime(datetime_str, '%Y-%m-%d at %H:%M:%S')
            .replace(tzinfo=datetime.timezone.utc)
            .timestamp()
        )

    def to_basis_set(self, basis_used: str):
        # checks necessary as LOBSTER 5.1.1 writes basis names now in lower case
        basis_used_lower = basis_used.lower()
        if basis_used_lower == 'pbevaspfit2015':
            positions = [3, 7]
        elif basis_used_lower in ['bunge', 'koga']:
            positions = [0]
        char_list = list(basis_used)
        # Loop through the positions and capitalize them if within bounds
        for pos in positions:
            if 0 <= pos < len(char_list):
                char_list[pos] = char_list[pos].upper()
        # Join the list back into a string
        return ''.join(char_list)

    def get_spilling(self, source: list[dict[str, Any]], type: str = '') -> np.ndarray:
        return np.array(
            [s.get(f'abs_{type}_spilling') for s in source], dtype=np.float64
        )


class LobsterStructureParser(MappingParser):
    code_name: str = ''

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self) -> ase.Atoms:
        code_name = self.code_name.lower()
        basedir = os.path.dirname(self.filepath)
        if code_name == 'vasp':
            structure_format = 'vasp'
            files = search_files('CONTCAR*', basedir)
        elif code_name == 'quantum espresso':
            structure_format = 'espresso-in'
            files = search_files('*scf.in', basedir)
        else:
            files = []

        if not files:
            raise self.logger.error('No Lobster structure file found.')

        return ase.io.read(files[0], format=structure_format)

    def to_dict(self):
        if self.data_object is not None:
            return self.data_object.todict()
        return {}

    def from_dict(self, data_dict: dict[str, Any]):
        pass

    def get_atoms(self):
        return [
            dict(number=number, symbol=self.data_object.symbols[n])
            for n, number in enumerate(self.data.get('numbers', []))
        ]


class LobsterTextParser(TextParser):
    _sources = []

    def to_dict(self) -> dict[str, Any]:
        dct = {}
        if self.data_object is None:
            return dct

        basedir = os.path.dirname(self.filepath)
        for source in self._sources:
            files = search_files(f'{source}.lobster*', basedir)
            if not files:
                continue
            self.data_object.mainfile = files[0]
            self.data_object.parse()
            dct[source] = self.data_object._results
        return dct


class LobsterICOXPLISTParser(LobsterTextParser):
    text_parser = ICOXPLISTParser()
    _sources = ['ICOHPLIST', 'ICOOPLIST', 'ICOBILIST']
    _atom_indices = []

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_dict(self):
        self._atom_indices = []
        return super().to_dict()

    def get_atom_pair_indices(self, source: dict[str, Any]) -> list[int]:
        if not self._atom_indices:
            data = source['data'][0]
            unique = []
            for n, index in enumerate(data[0]):
                if index not in unique:
                    unique.append(index)
                    self._atom_indices.append(n)
        return self._atom_indices

    def get_orbital_pairs(self, index: int, source: np.ndarray) -> list[int]:
        data = source['data'][0]
        atom_index = data[0][index]
        return [
            dict(atomMU=data[1][n], atomNU=data[2][n])
            for n, a_n in enumerate(data[0])
            if a_n == atom_index
        ][1:]

    def get_bond_label(self, n: int, source: np.ndarray) -> str:
        return f'{source[1][n]} -> {source[2][n]}'

    def get_atom_labels(self, source: dict[str, Any], atom: int = 0) -> list[str]:
        labels = [
            source['data'][0][atom + 1][n] for n in self.get_atom_pair_indices(source)
        ]
        return labels

    def get_atom_distances(self, source: dict[str, Any]) -> list[float]:
        return np.array(
            [source['data'][0][3][n] for n in self.get_atom_pair_indices(source)],
            dtype=np.float64,
        )

    def get_translations(self, source: dict[str, Any]) -> np.ndarray | None:
        data = source['data'][0]
        if data[4][0].lstrip('-').isdigit():
            # return np.array(data[4:7], dtype=np.int32).T
            data = data.T
            indices = self.get_atom_pair_indices(source)
            return np.array(data[indices][:, 4:7], dtype=np.int32)
        return None

    def get_integrated_coxp_at_fermi_level(self, source: dict[str, Any]) -> np.ndarray:
        icoxp = []
        indices = self.get_atom_pair_indices(source)
        line_split = source['data'][0].T[0]
        case = 0
        if len(line_split) == 9 and not line_split[-1].isdigit():  # noqa: PLR2004
            # Spin polarized data LOBSTER version 5.1 and above
            case = 1
        elif len(line_split) == 6 and line_split[-1].isdigit():  # noqa: PLR2004
            case = 2
        for data in source['data']:
            if case == 1:
                icoxp.append(data[7][indices])
                icoxp.append(data[8][indices])
            elif case == 2:  # noqa: PLR2004
                icoxp.append(data[4][indices])
            else:
                icoxp.append(data[-1][indices])
        return np.array(icoxp, dtype=np.float64)

    def get_bonds(self, source: dict[str, Any]) -> np.ndarray | None:
        data = source['data'][0]
        if data[-1][0].isdigit():
            return np.array(data[-1], dtype=np.int32)
        return None


class LobsterCOXPCARParser(LobsterTextParser):
    _sources = ['COHPCAR', 'COOPCAR', 'COBICAR']
    _nspin = 0

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self) -> FileParser:
        return COXPCARParser()

    def to_dict(self):
        self._nspin = 0
        return super().to_dict()

    def get_bond_label(self, source: tuple) -> str:
        return f'{source[1]} -> {source[2]}'

    def get_orbital_pairs(
        self, bond_pair: tuple, source: dict[str, Any], name=''
    ) -> list[dict[str, Any]]:
        coxp_pairs = source.get('coxp_pairs', [])
        return [n for n, pair in enumerate(coxp_pairs) if pair[0] == bond_pair[0]][1:]

    def get_atom_label(self, n: int, source: dict[str, Any], atom: int = 0) -> str:
        return (
            source.get('coxp_pairs', [])[n][atom + 1].replace('[', '_').replace(']', '')
        )

    def get_atom_value(self, source: np.ndarray, coxp_pairs: list[tuple]) -> np.ndarray:
        indices = [
            n
            for n, pair in enumerate(coxp_pairs)
            if '[' not in pair[1] and '[' not in pair[2]
        ]
        return source[indices]

    def get_orbital_value(
        self, n: int, source: dict[str, Any], type: str = 'pair_icoxp', name=''
    ) -> np.ndarray | None:
        return source.get(type, [])[n]


class LobsterCHARGEParser(LobsterTextParser):
    text_parser = CHARGEParser()
    _sources = ['CHARGE']

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self) -> TextParser:
        charge_files = search_files('CHARGE.lobster*', os.path.dirname(self.filepath))
        if not charge_files:
            self.logger.error('No Lobster CHARGE file found.')
            return self.text_parser
        self.text_parser.mainfile = charge_files[0]
        return self.text_parser

    def get_contributions(
        self, source: dict[str, Any], kind: str
    ) -> list[dict[str, Any]]:
        contributions = []
        charges = source.get(kind)
        for n, symbol in enumerate(source.get('symbols', [])):
            contribution = dict(symbol=symbol)
            if charges is not None:
                contribution['value'] = charges[n]
            contributions.append(contribution)
        return contributions

    def get_charges(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        charges = []
        total = source.get('total')
        for n, kind in enumerate(self.text_parser.kinds):
            charge = dict(kind=kind, value=total[n])
            contributions = self.get_contributions(source.get('charges', {}), kind)
            if contributions:
                charge['contributions'] = contributions
            charges.append(charge)
        return charges


class LobsterDOSCARParser(LobsterTextParser):
    _sources = ['DOSCAR']

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self) -> FileParser:
        return DOSCARParser()

    def to_dict(self):
        doscar_files = search_files(
            'DOSCAR.LSO.lobster*', os.path.dirname(self.filepath)
        )
        if doscar_files:
            self._sources = ['DOSCAR.LSO']
        return super().to_dict()

    def get_dos(
        self, total: np.ndarray, projected: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        dos = []
        n_spin = len(total) // 2
        for spin in range(n_spin):
            dct = dict(
                dos=total[spin],
                integrated=total[n_spin + spin],
                spin=spin,
                projected=[],
            )
            for atom, pdos_dict in enumerate(projected):
                for lm, pdos_lm in enumerate(pdos_dict['dos'][spin]):
                    dct['projected'].append(
                        dict(dos=pdos_lm, lm=lm, atom_index=atom, spin=spin)
                    )
            dos.append(dct)
        return dos


class LobsterMetainfoParser(MetainfoParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class LobsterArchiveWriter(ArchiveWriter):
    code_name = 'LOBSTER'
    mainfile_parser = LobsterMainfileParser()
    metainfo_parser = LobsterMetainfoParser()
    structure_parser = LobsterStructureParser()
    icoxplist_parser = LobsterICOXPLISTParser()
    coxpcar_parser = LobsterCOXPCARParser()
    charge_parser = LobsterCHARGEParser()
    doscar_parser = LobsterDOSCARParser()

    def write_to_archive(self):
        self.archive.data = Simulation(program=Program(name=self.code_name))
        self.metainfo_parser.data_object = self.archive.data

        # parser mainfile
        self.metainfo_parser.annotation_key = lobster.OUT_KEY
        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.convert(self.metainfo_parser)

        # parse structure
        self.metainfo_parser.annotation_key = lobster.STRUCTURE_KEY
        self.structure_parser.code_name = self.mainfile_parser.data.get(
            'x_lobster_code'
        )
        self.structure_parser.filepath = self.mainfile
        self.structure_parser.convert(self.metainfo_parser)

        # parse ICOXPLIST
        self.metainfo_parser.annotation_key = lobster.ICOXPLIST_KEY
        self.icoxplist_parser.filepath = self.mainfile
        self.icoxplist_parser.convert(self.metainfo_parser)

        # parse COXPCAR
        self.metainfo_parser.annotation_key = lobster.COXPCAR_KEY
        self.coxpcar_parser.filepath = self.mainfile
        self.coxpcar_parser.convert(self.metainfo_parser)

        # parse CHARGE
        self.metainfo_parser.annotation_key = lobster.CHARGE_KEY
        self.charge_parser.filepath = self.mainfile
        self.charge_parser.convert(self.metainfo_parser)

        # parse DOSCAR
        self.metainfo_parser.annotation_key = lobster.DOSCAR_KEY
        self.doscar_parser.filepath = self.mainfile
        self.doscar_parser.convert(self.metainfo_parser)


class LobsterParser(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger) -> None:
        archive_writer = LobsterArchiveWriter()
        archive_writer.write(mainfile, archive, logger)
