import datetime
import os
import re
from importlib import reload
from typing import Any

import numpy as np
import pint
from ase.data import chemical_symbols
from nomad import atomutils
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import crystal

from .file_parser import F25Parser, OutputParser


class CrystalOutputParser(TextParser):
    libxc_map = {
        'PBEXC': ['GGA_C_PBE', 'GGA_X_PBE'],
        'PBE0': ['HYB_GGA_XC_PBEH'],
        'B3LYP': ['HYB_GGA_XC_B3LYP'],
        'HSE06': ['HYB_GGA_XC_HSE06'],
        'M06': ['MGGA_C_M06', 'HYB_MGGA_X_M06'],
        'M05-2X': ['HYB_MGGA_XC_M05_2X'],
        'LC-WPBE': ['HYB_GGA_XC_LRC_WPBE'],
        'PBE': ['GGA_X_PBE', 'GGA_C_PBE'],
        'PBESOL': ['GGA_X_PBE_SOL', 'GGA_C_PBE_SOL'],
        'BECKE': ['GGA_X_B88'],
        'LDA': ['LDA_X'],
        'PWGGA': ['GGA_X_PW91', 'GGA_C_PW91'],
        'PZ': ['LDA_C_PZ'],
        'WFN': ['LDA_C_VWN'],
    }

    @property
    def logger(self):
        pass

    def to_unix_time(self, value: str) -> float | None:
        """Transforms the Crystal-specific float notation into a floating point
        number.
        """
        if value is None:
            return None

        value = value.strip()
        date_time_obj = datetime.datetime.strptime(value, '%d %m %Y TIME %H:%M:%S.%f')
        return date_time_obj.timestamp()

    def get_lattice_vectors(self, source: dict[str, Any]) -> pint.Quantity:
        lattice_vectors = source.get(
            'lattice_parameters', source.get('lattice_vectors_restart')
        )
        if lattice_vectors.shape == (6,):
            lattice_vectors = atomutils.cellpar_to_cell(lattice_vectors, degrees=True)
        return lattice_vectors * ureg.angstrom

    def get_positions(self, source: dict[str, Any]) -> pint.Quantity:
        labels_positions = source.get(
            'labels_positions_nanotube', source.get('labels_positions')
        )

        if labels_positions is None:
            labels_positions = self.data.get('labels_positions_restart')

        if labels_positions is None:
            return labels_positions

        positions = labels_positions[:, 4:7].astype(np.float64)

        dimensionality = self.data.get('dimensionality')
        lattice_vectors = self.get_lattice_vectors(source)

        scaled_pos = np.zeros((len(positions), 3), dtype=np.float64)
        scaled_pos[:, :dimensionality] = positions[:, :dimensionality]
        if lattice_vectors is not None:
            cart_pos = atomutils.to_cartesian(scaled_pos, lattice_vectors.magnitude)
            cart_pos[:, dimensionality:] = positions[:, dimensionality:]
        else:
            cart_pos = scaled_pos

        return cart_pos * ureg.angstrom

    def get_atoms(self, source: dict[str, Any]) -> list[str]:
        labels_positions = source.get(
            'labels_positions_nanotube', source.get('labels_positions')
        )
        labels = []
        numbers = []
        if labels_positions is not None:
            labels = labels_positions[:, 3]
            numbers = labels_positions[:, 2]
        else:
            labels_positions = self.data.get('labels_positions_restart')
            if labels_positions is not None:
                labels = labels_positions[:, 2]
                numbers = labels_positions[:, 1]

        def normalize_label(label: str) -> str:
            norm = re.match(r'([a-z][a-z]?).*', label.lower())
            # unknown specie
            # TODO not possible to define ghost atom
            unknown = None
            if norm:
                label = norm.group(1).capitalize()
                return label if label in chemical_symbols[1:] else unknown
            return unknown

        return [
            dict(label=normalize_label(label), number=numbers[n])
            for n, label in enumerate(labels)
        ] or None

    def get_xc_functionals(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        xc_functionals = set()
        for key in ['exchange', 'correlation', 'exchange_correlation']:
            xc_functionals.update(self.libxc_map.get(source.get(key), []))
        return [dict(name=xc) for xc in xc_functionals]

    def get_outputs(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        outputs = []
        for n, step in enumerate(source.get('geo_opt', {}).get('geo_opt_step', [])):
            outputs.append(dict(energy=step.get('energy')))

        output0 = outputs[0] if outputs else dict()
        output0.setdefault('energy', source.get('energy_total'))
        forces = source.get('forces')
        if forces is not None:
            output0.setdefault(
                'forces', forces[:, 2:].astype(float) * ureg.hartree / ureg.bohr
            )

        return outputs

    def get_systems(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        initial = source.get('system_edited', source)
        systems = [
            # initial system
            dict(
                positions=self.get_positions(initial),
                atoms=self.get_atoms(initial),
                lattice_vectors=self.get_lattice_vectors(initial),
            )
        ]
        # skip first step same as initial
        for step in source.get('geo_opt', {}).get('geo_opt_step', [])[1:]:
            systems.append(
                dict(
                    positions=self.get_positions(step),
                    atoms=self.get_atoms(step),
                    lattice_vectors=self.get_lattice_vectors(step),
                )
            )
        return systems

    def get_band_structures(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        pass


class CrystalF25Parser(TextParser):
    @property
    def logger(self):
        pass

    @staticmethod
    def to_array(cols: int, rows: int, values: str) -> np.ndarray:
        """Transforms the Crystal-specific f25 array syntax into a numpy array."""
        values = values.replace('\n', '').replace('\r', '')
        values = [values[n : n + 12] for n in range(0, len(values), 12)]
        # does not seem to work
        # values = textwrap.wrap(values, 12)
        return np.array(values, dtype=np.float64).reshape((rows, cols))

    def get_dos(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        first_row = source['first_row']
        cols, rows = [int(first_row[n]) for n in range(2)]
        de = first_row[3]
        # fermi_energy = first_row[4]
        second_row = source['second_row']
        start_energy = second_row[1]
        dos_values = self.to_array(cols, rows, source['values']).T
        return [
            dict(
                energies=(start_energy + np.arange(rows) * de) * ureg.hartree,
                values=dos_values[n],
            )
            for n in range(len(dos_values))
        ]


class CrystalMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        pass


class CrystalArchiveWriter(ArchiveWriter):
    output_parser = CrystalOutputParser(text_parser=OutputParser())
    f25_parser = CrystalF25Parser(text_parser=F25Parser())
    archive_parser = CrystalMetainfoParser()

    def write_to_archive(self):
        # reload schema to update annotations
        reload(crystal)

        # main output file
        self.archive_parser.annotation_key = 'out'
        self.archive_parser.data_object = Simulation(program=Program(name='Crystal'))

        self.output_parser.filepath = self.mainfile
        self.output_parser.convert(self.archive_parser)

        self.archive.data = self.archive_parser.data_object

        f25_filepath = self.output_parser.data.get(
            'f25_filepath1', self.output_parser.data.get('f25_filepath2')
        )
        if f25_filepath:
            # parser f25 file
            self.archive_parser.annotation_key = 'f25'

            self.f25_parser.filepath = os.path.join(
                os.path.dirname(self.mainfile), os.path.basename(f25_filepath)
            )

            self.f25_parser.convert(self.archive_parser)


class CrystalParser(MatchingParser):
    archive_writer = CrystalArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_arhives: dict[str, EntryArchive] = {},
    ):
        self.archive_writer.write(mainfile, archive, logger, child_arhives)
