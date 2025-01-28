import os
import re
from collections.abc import Iterable
from typing import Any, Union

import numpy as np
from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import TaskReference
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import (
    MetainfoParser,
)
from nomad.parsing.file_parser.mapping_parser import (
    TextParser as TextMappingParser,
)
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import DFTGWWorkflow, SinglePoint
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.fhiaims.out_parser import (
    RE_GW_FLAG,
    FHIAimsOutFileParser,
)
from nomad_simulation_parsers.parsers.utils.general import remove_mapping_annotations


class FHIAimsOutMappingParser(TextMappingParser):
    _gw_flag_map = {
        'gw': 'G0W0',
        'gw_expt': 'G0W0',
        'ev_scgw0': 'ev-scGW',
        'ev_scgw': 'ev-scGW',
        'scgw': 'scGW',
    }

    _xc_map = {
        'Perdew-Wang parametrisation of Ceperley-Alder LDA': [
            {'name': 'LDA_C_PW'},
            {'name': 'LDA_X'},
        ],
        'Perdew-Zunger parametrisation of Ceperley-Alder LDA': [
            {'name': 'LDA_C_PZ'},
            {'name': 'LDA_X'},
        ],
        'VWN-LDA parametrisation of VWN5 form': [
            {'name': 'LDA_C_VWN'},
            {'name': 'LDA_X'},
        ],
        'VWN-LDA parametrisation of VWN-RPA form': [
            {'name': 'LDA_C_VWN_RPA'},
            {'name': 'LDA_X'},
        ],
        'AM05 gradient-corrected functionals': [
            {'name': 'GGA_C_AM05'},
            {'name': 'GGA_X_AM05'},
        ],
        'BLYP functional': [{'name': 'GGA_C_LYP'}, {'name': 'GGA_X_B88'}],
        'PBE gradient-corrected functionals': [
            {'name': 'GGA_C_PBE'},
            {'name': 'GGA_X_PBE'},
        ],
        'PBEint gradient-corrected functional': [
            {'name': 'GGA_C_PBEINT'},
            {'name': 'GGA_X_PBEINT'},
        ],
        'PBEsol gradient-corrected functionals': [
            {'name': 'GGA_C_PBE_SOL'},
            {'name': 'GGA_X_PBE_SOL'},
        ],
        'RPBE gradient-corrected functionals': [
            {'name': 'GGA_C_PBE'},
            {'name': 'GGA_X_RPBE'},
        ],
        'revPBE gradient-corrected functionals': [
            {'name': 'GGA_C_PBE'},
            {'name': 'GGA_X_PBE_R'},
        ],
        'PW91 gradient-corrected functionals': [
            {'name': 'GGA_C_PW91'},
            {'name': 'GGA_X_PW91'},
        ],
        'M06-L gradient-corrected functionals': [
            {'name': 'MGGA_C_M06_L'},
            {'name': 'MGGA_X_M06_L'},
        ],
        'M11-L gradient-corrected functionals': [
            {'name': 'MGGA_C_M11_L'},
            {'name': 'MGGA_X_M11_L'},
        ],
        'TPSS gradient-corrected functionals': [
            {'name': 'MGGA_C_TPSS'},
            {'name': 'MGGA_X_TPSS'},
        ],
        'TPSSloc gradient-corrected functionals': [
            {'name': 'MGGA_C_TPSSLOC'},
            {'name': 'MGGA_X_TPSS'},
        ],
        'hybrid B3LYP functional': [{'name': 'HYB_GGA_XC_B3LYP5'}],
        'Hartree-Fock': [{'name': 'HF_X'}],
        'HSE': [{'name': 'HYB_GGA_XC_HSE03'}],
        'HSE-functional': [{'name': 'HYB_GGA_XC_HSE06'}],
        'hybrid-PBE0 functionals': [
            {'name': 'GGA_C_PBE'},
            {
                'name': 'GGA_X_PBE',
                'weight': lambda x: 0.75 if x is None else 1.0 - x,
            },
            {'name': 'HF_X', 'weight': lambda x: 0.25 if x is None else x},
        ],
        'hybrid-PBEsol0 functionals': [
            {'name': 'GGA_C_PBE_SOL'},
            {
                'name': 'GGA_X_PBE_SOL',
                'weight': lambda x: 0.75 if x is None else 1.0 - x,
            },
            {'name': 'HF_X', 'weight': lambda x: 0.25 if x is None else x},
        ],
        'Hybrid M06 gradient-corrected functionals': [
            {'name': 'MGGA_C_M06'},
            {'name': 'HYB_MGGA_X_M06'},
        ],
        'Hybrid M06-2X gradient-corrected functionals': [
            {'name': 'MGGA_C_M06_2X'},
            {'name': 'HYB_MGGA_X_M06'},
        ],
        'Hybrid M06-HF gradient-corrected functionals': [
            {'name': 'MGGA_C_M06_HF'},
            {'name': 'HYB_MGGA_X_M06'},
        ],
        'Hybrid M08-HX gradient-corrected functionals': [
            {'name': 'MGGA_C_M08_HX'},
            {'name': 'HYB_MGGA_X_M08_HX'},
        ],
        'Hybrid M08-SO gradient-corrected functionals': [
            {'name': 'MGGA_C_M08_SO'},
            {'name': 'HYB_MGGA_X_M08_SO'},
        ],
        'Hybrid M11 gradient-corrected functionals': [
            {'name': 'MGGA_C_M11'},
            {'name': 'HYB_MGGA_X_M11'},
        ],
    }

    _section_names = ['full_scf', 'geometry_optimization', 'molecular_dynamics']

    def get_fhiaims_file(self, default: str) -> list[str]:
        maindir = os.path.dirname(self.filepath)
        base, *ext = default.split('.')
        ext = '.'.join(ext)
        base = base.lower()
        files = os.listdir(maindir)
        files = [os.path.basename(f) for f in files]
        files = [
            os.path.join(maindir, f)
            for f in files
            if base.lower() in f.lower() and f.endswith(ext)
        ]
        files.sort()
        return files

    def get_xc_functionals(self, xc: str) -> list[dict[str, Any]]:
        return [
            dict(name=functional.get('name')) for functional in self._xc_map.get(xc, [])
        ]

    def get_dos(
        self,
        total_dos_files: list[list[str]],
        atom_dos_files: list[list[str]],
        species_dos_files: list[list[str]],
    ) -> list[dict[str, Any]]:
        def load_dos(dos_file: str) -> list[dict[str, Any]]:
            dos_files = self.get_fhiaims_file(dos_file)
            if not dos_files:
                return []
            try:
                data = np.loadtxt(dos_files[0]).T
            except Exception:
                return []
            if not np.size(data):
                return []
            return [
                dict(energies=data[0], values=value, nenergies=len(data[0]))
                for value in data[1:]
            ]

        def get_pdos(dos_files: list[str], dos_labels: list[str], type=str):
            dos = []
            for dos_file in dos_files:
                labels = [label for label in dos_labels if label in dos_file]
                pdos = load_dos(dos_file)
                if not pdos:
                    continue
                for n, data in enumerate(pdos):
                    # TODO use these to link pdos to system
                    data['type'] = type
                    data['label'] = labels[n % len(labels)]
                    data['spin'] = 1 if 'spin_dn' in dos_file else 0
                    data['orbital'] = n - 1 if n else None
                dos.extend(pdos)
            return dos

        projected_dos = []
        # atom-projected dos
        if atom_dos_files:
            projected_dos.extend(get_pdos(*atom_dos_files, type='atom'))

        # species-projected dos
        if species_dos_files:
            projected_dos.extend(get_pdos(*species_dos_files, type='species'))

        # total dos
        total_dos = []
        for dos_file in (
            total_dos_files[0] if total_dos_files else ['KS_DOS_total_raw.dat']
        ):
            dos = load_dos(dos_file)
            for n, data in enumerate(dos):
                data['spin'] = n
                pdata = data.setdefault('projected_dos', [])
                pdata.extend([d for d in projected_dos if d['spin'] == data['spin']])
            total_dos.extend(dos)

        return total_dos

    def get_eigenvalues(
        self, source: list[dict[str, Any]], params: dict[str, Any]
    ) -> list[dict[str, Any]]:
        n_spin = params.get('Number of spin channels', 1)
        eigenvalues = []
        for data in source:
            kpts = data.get('kpoints', [np.zeros(3)] * n_spin)
            kpts = np.reshape(kpts, (len(kpts) // n_spin, n_spin, 3))
            kpts = np.transpose(kpts, axes=(1, 0, 2))[0]

            occs_eigs = data.get('occupation_eigenvalue')
            n_kpts = len(kpts)
            n_eigs = len(occs_eigs) // (n_kpts * n_spin)
            occs_eigs = np.transpose(
                np.reshape(occs_eigs, (n_kpts, n_spin, n_eigs, 2)), axes=(3, 1, 0, 2)
            )
            for spin in range(n_spin):
                eigenvalues.append(
                    dict(
                        nbands=n_eigs,
                        npoints=n_kpts,
                        points=kpts,
                        occupations=occs_eigs[0][spin],
                        eigenvalues=occs_eigs[1][spin],
                    )
                )
        return eigenvalues

    def get_energies(self, source: dict[str, Any]) -> dict[str, Any]:
        total_keys = ['Total energy uncorrected', 'Total energy']
        energies = {}
        components = []
        for key, val in source.get('energy', {}).items():
            if key in total_keys:
                energies.setdefault('value', val)
            else:
                components.append({'name': key, 'value': val})
        for key, val in source.get('energy_components', [{}])[-1].items():
            components.append({'name': key, 'value': val})
        energies['components'] = components
        return energies

    def get_forces(self, source: dict[str, Any]) -> dict[str, Any]:
        return dict(
            forces=source.get('forces'), npoints=len(source.get('forces', [])), rank=[3]
        )

    def get_gw_flag(self, gw_flag: str):
        return self._gw_flag_map.get(gw_flag)

    def get_sections(self, source: dict[str, Any], **kwargs) -> list[dict[str, Any]]:
        result = []
        include = kwargs.get('include')
        for name in self._section_names:
            for data in source.get(name, []):
                res = {}
                for key in data.keys():
                    if include and key not in include:
                        continue
                    val = data.get(key, self.data.get(key))
                    if val is not None:
                        res[key] = val
                if res:
                    result.append(res)
        return result


class FHIAimsArchiveWriter(ArchiveWriter):
    annotation_key: str = 'text'

    def write_to_archive(
        self,
    ) -> None:
        from nomad_simulation_parsers.schema_packages import fhiaims

        out_parser = FHIAimsOutMappingParser()
        out_parser.text_parser = FHIAimsOutFileParser()
        out_parser.filepath = self.mainfile

        archive_data_handler = MetainfoParser()
        archive_data_handler.annotation_key = self.annotation_key
        archive_data_handler.data_object = Simulation(program=Program(name='FHI-aims'))

        out_parser.convert(archive_data_handler, remove=True)

        self.archive.data = archive_data_handler.data_object

        # separate parsing of dos due to a problem with mapping physical
        # property variables
        archive_data_handler.annotation_key = 'text_dos'
        out_parser.convert(archive_data_handler, remove=True)

        gw_archive = self.child_archives.get('GW') if self.child_archives else None
        if gw_archive is not None:
            self.archive.workflow2 = SinglePoint(name='DFT')

            # GW single point
            parser = FHIAimsArchiveWriter()
            parser.annotation_key = 'text_gw'
            parser.write(self.mainfile, gw_archive, self.logger)
            gw_archive.workflow2 = SinglePoint(name='GW')

            # DFT-GW workflow
            gw_workflow_archive = self.child_archives.get('GW_workflow')
            gw_workflow_archive.workflow2 = DFTGWWorkflow(tasks=[
                TaskReference(task=self.archive.workflow2),
                TaskReference(task=gw_archive.workflow2)
            ])

        # close file contexts
        self.out_parser = out_parser
        # out_parser.close()
        archive_data_handler.close()

        # remove annotations
        remove_mapping_annotations(fhiaims.general.Simulation.m_def)


class FHIAimsParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = FHIAimsArchiveWriter()

    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> Union[bool, Iterable[str]]:
        is_mainfile = super().is_mainfile(
            filename=filename,
            mime=mime,
            buffer=buffer,
            decoded_buffer=decoded_buffer,
            compression=compression,
        )
        if is_mainfile:
            match = re.search(RE_GW_FLAG, decoded_buffer)
            if match:
                gw_flag = match[1]
            else:
                gw_flag = None
                with open(filename) as f:
                    while True:
                        line = f.readline()
                        match = re.match(RE_GW_FLAG, f'\n{line}')
                        if match:
                            gw_flag = match[1]
                            break
                        if not line:
                            break
            if gw_flag in FHIAimsOutMappingParser._gw_flag_map.keys():
                self.creates_children = True
                return ['GW', 'GW_workflow']
        return is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, EntryArchive] = None,
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
