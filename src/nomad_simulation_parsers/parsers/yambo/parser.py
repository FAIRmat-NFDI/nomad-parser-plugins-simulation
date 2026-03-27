import os
from datetime import datetime
from importlib import reload
from typing import Any

import numpy as np
from ase.data import chemical_symbols
from netCDF4 import Dataset
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import (
    MappingParser,
    MetainfoParser,
    TextParser,
)
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import yambo

from .file_parsers import MainfileParser, NetCDFParser

LOGGER = get_logger(__name__)


# TODO temporary fix for structlog unable to propagate logger
class YamboMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class YamboNetCDFParser(MappingParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.netcdf_parser = NetCDFParser(mainfile=self.filepath)

    @property
    def logger(self):
        return LOGGER

    def load_file(self) -> Dataset:
        self.netcdf_parser.mainfile = self.filepath
        return self.netcdf_parser.netcdf_file

    def to_dict(self, **kwargs) -> dict[str, Any]:
        self.netcdf_parser.parse()
        return self.netcdf_parser.results

    def from_dict(self, dct: dict[str, Any]) -> None:
        pass

    def get_positions(self) -> np.ndarray | None:
        positions = self.data.get('ATOM_POS', [])
        max_n_atoms = self.data.get('MAX_ATOMS', [0])[0]
        n_atoms = self.data.get('N_ATOMS', [])
        if not max_n_atoms or not n_atoms or len(positions) == 0:
            return None
        # We split the positions array into blocks, each corresponding
        # to a chemical species, we extract the first n_atoms only
        # (value of n_atoms for each chemical species present in the system)
        # from each block, and we reassemble the modified blocks
        # into the corrected positions array
        positions = np.array(positions)
        selected = []
        positions = positions.reshape(-1, 3)
        n_points = positions.shape[0]
        n_blocks = int(n_points / max_n_atoms)

        if n_blocks != len(n_atoms):
            return None

        for i in range(n_blocks):
            start_idx = int(i * max_n_atoms)
            end_idx = int((i + 1) * max_n_atoms)
            block = positions[start_idx:end_idx]
            n_to_select = int(n_atoms[int(i)])
            selected_from_block = block[:n_to_select]
            for point in selected_from_block:
                selected.append(point)
        return np.array(selected)

    def get_labels(self) -> list[str]:
        n_atoms = self.data.get('N_ATOMS', [])
        atomic_numbers = self.data.get('atomic_numbers', [])
        if not n_atoms or not atomic_numbers:
            return []
        atom_numbers = np.hstack(
            [
                [atomic_numbers[int(n)]] * int(n_atoms[int(n)])
                for n in range(len(n_atoms))
            ]
        )
        return [chemical_symbols[int(num)] for num in atom_numbers]

    def get_lattice_vectors(self) -> np.ndarray:
        lattice_vectors = np.array(self.data.get('LATTICE_VECTORS'), dtype=np.float64)
        return lattice_vectors

    def get_kpoints(self) -> np.ndarray | None:
        if self.data.get('K-POINTS') is not None:
            return self.data.get('K-POINTS').T
        elif self.data.get('QP_kpts') is not None:
            return self.data.get('QP_kpts').T
        else:
            return None

    def get_eigenvalues(self) -> list[dict[str, np.ndarray]]:
        eigenvalues = []
        qp_table = self.data.get('QP_table')
        n_spin = qp_table.shape[1] // 2 if qp_table is not None else 0
        if (qp_e := self.data.get('EIGENVALUES')) is not None:
            eigenvalues.extend(
                [dict(energies=eig * ureg.eV) for n, eig in enumerate(qp_e)]
            )
        if (
            (qp_e_eo_z := self.data.get('QP_E_Eo_Z')) is not None
            or (qp_e := self.data.get('QP_E')) is not None
            and n_spin
        ):
            if qp_e_eo_z is not None:
                qp_energy, bare_energy, z = qp_e_eo_z[0].T
            else:
                qp_energy = qp_e.T[0]
                bare_energy = self.data.get('QP_Eo')
                z = self.data.get('QP_Z').T[0]
            # TODO verify if indeed energies are only for one kpoint
            shape = (n_spin, 1, len(qp_energy) // n_spin)
            value_qp = np.reshape(qp_energy, shape) * ureg.hartree
            value_ks = np.reshape(bare_energy, shape) * ureg.hartree
            qp_linearization_prefactor = np.reshape(z, shape)
            eigenvalues.extend(
                [
                    dict(
                        value_qp=value_qp[n],
                        value_ks=value_ks[n],
                        qp_linearization_prefactor=qp_linearization_prefactor[n],
                    )
                    for n in range(n_spin)
                ]
            )
        if (
            (sx_vxc := self.data.get('Sx_Vxc')) is not None
            or (sx := self.data.get('Sx')) is not None
            and n_spin
        ):
            if sx_vxc is not None:
                if sx_vxc.shape[0] % 8 == 0:
                    qp = sx_vxc.reshape(-1, 8).T
                    sx, vxc = qp[4], qp[6]
                else:
                    qp = sx_vxc.reshape(-1, 7).T
                    sx, vxc = qp[3], qp[5]
            else:
                sx = sx.T[0]
                vxc = self.data.get('Vxc').T[0]
            shape = (n_spin, 1, len(sx) // n_spin)
            value_exchange = np.reshape(sx, shape)
            value_xc_potential = np.reshape(vxc, shape)
            eigenvalues.extend(
                [
                    dict(
                        value_exchange=value_exchange[n] * ureg.hartree,
                        value_xc_potential=value_xc_potential[n] * ureg.hartree,
                    )
                    for n in range(n_spin)
                ]
            )
        return eigenvalues


class YamboMainfileParser(TextParser):
    @property
    def logger(self):
        return LOGGER

    def get_wallstart(self, parsed: str) -> float:
        return datetime.strptime(parsed, '%d/%m/%Y %H:%M').timestamp()

    def get_outputs(
        self,
        energies_occupations: dict[str, Any],
        modules: list[dict[str, Any]],
        transferred_momenta: dict[str, Any],
    ) -> list[dict[str, Any]]:
        outputs = []
        data = {}
        for key, val in energies_occupations.items():
            if key == 'eigenenergies':
                kpoints = val.get('kpoints')
                energies = val.get('energies')
                if kpoints is None or energies is None:
                    continue
                # TODO deal with spin polarized data
                n_spin = 1
                energies = (
                    np.reshape(
                        energies,
                        (
                            n_spin,
                            len(kpoints),
                            np.size(energies) // len(kpoints),
                        ),
                    )
                    * ureg.eV
                )
                data['eigenvalues'] = [
                    dict(energies=energies[n]) for n in range(n_spin)
                ]
            else:
                data[key] = val
        if data:
            outputs.append(data)

        def get_qp_properties(source: dict[str, Any]) -> list[dict[str, Any]]:
            qp_energy = source.get('qp_energy')
            if qp_energy is None:
                return []
            energies = np.transpose([q.band for q in qp_energy])
            qp_energy = energies[2].T
            n_spin = 1
            value_qp = np.reshape(qp_energy, (n_spin, *np.shape(qp_energy))) * ureg.eV
            value_ks = (
                np.reshape(energies[1].T, (n_spin, *np.shape(qp_energy))) * ureg.eV
            )
            qp_linearization_prefactor = np.reshape(
                energies[4].T, (n_spin, *np.shape(qp_energy))
            )
            return [
                dict(
                    value_qp=value_qp[n],
                    value_ks=value_ks[n],
                    qp_linearization_prefactor=qp_linearization_prefactor[n],
                )
                for n in range(n_spin)
            ]

        def unpack_modules(modules: list[dict[str, Any]]) -> list[dict[str, Any]]:
            module_names = ['dyson', 'local_xc_nonlocal_fock', 'bare_xc']
            outputs = []
            for module in modules:
                for name in module_names:
                    source = module.get(name, {})
                    data = {}
                    for key, val in source.items():
                        if key == 'qp_properties':
                            data['eigenvalues'] = get_qp_properties(val)
                        else:
                            data[key] = val
                    if data:
                        outputs.append(data)
            return outputs

        outputs.extend(unpack_modules(modules or []))

        if qp_properties := get_qp_properties(
            transferred_momenta.get('qp_properties', {})
        ):
            outputs.append(dict(eigenvalues=qp_properties))

        outputs.extend(unpack_modules(transferred_momenta.get('modules', [])))

        return outputs

    def get_band_gaps(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        if not hasattr(source, 'get'):
            return []

        required_levels = 2
        valence_conduction = source.get('valence_conduction')
        has_valence_conduction = (
            valence_conduction is not None
            and len(valence_conduction) >= required_levels
        )
        if has_valence_conduction:
            valence, conduction = valence_conduction[0], valence_conduction[1]
        else:
            valence = source.get('valence')
            conduction = source.get('conduction')

        if valence is None or conduction is None:
            return []

        value = conduction - valence
        if hasattr(value, 'magnitude'):
            if value.magnitude < 0:
                value = 0 * value.units
        else:
            value = max(0.0, value)

        return [dict(value=value)]


class YamboArchiveWriter(ArchiveWriter):
    def write_to_archive(self):
        from nomad_simulation_parsers.schema_packages.yambo import Simulation

        data = Simulation()

        self.archive.data = data

        # set up parser for simulation data
        data_parser = YamboMetainfoParser()
        data_parser.data_object = data

        # set up parser for yambo main file
        mainfile_parser = YamboMainfileParser(text_parser=MainfileParser())
        mainfile_parser.filepath = self.mainfile

        # map mainfile data to simulation
        data_parser.annotation_key = yambo.OUT_KEY
        mainfile_parser.convert(data_parser)

        netcdf_file = (
            mainfile_parser.data.get('cpu_files_io', {})
            .get('input', {})
            .get('file', '')
        )
        netcdf_parser = None
        if netcdf_file:
            # set up parser for yambo netcdf file
            netcdf_parser = YamboNetCDFParser(
                filepath=os.path.join(os.path.dirname(self.mainfile), netcdf_file)
            )
            data_parser.annotation_key = yambo.NETCDF_KEY
            netcdf_parser.convert(data_parser)

        data_parser.close()
        if netcdf_parser is not None:
            netcdf_parser.close()


class YamboParser(MatchingParser):
    archive_writer = YamboArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        # reload schema to load yambo annotations
        reload(yambo)

        self.archive_writer.write(mainfile, archive, logger, child_archives)
