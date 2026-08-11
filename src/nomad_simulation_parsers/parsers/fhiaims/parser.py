import os
import re
from collections.abc import Iterable
from typing import Any

import numpy as np
from ase import Atoms
from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Link, TaskReference
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser
from nomad_file_parser.mapping_parser import TextParser as TextMappingParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import (
    DFTGWWorkflow,
    GeometryOptimization,
    MolecularDynamics,
    Phonon,
    SinglePoint,
)
from nomad_simulations.schema_packages.workflow.general import (
    EnergyConvergenceTarget,
    ForceConvergenceTarget,
    SimulationTaskReference,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.single_point import SinglePointMethod
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.fhiaims.out_parser import (
    RE_GW_FLAG,
    FHIAimsOutFileParser,
)
from nomad_simulation_parsers.parsers.phonopy.parser import phonopy_obj_to_archive
from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import fhiaims

from .common import ControlParser, GeometryParser

LOGGER = get_logger(__name__)


# TODO temporary fix for structlog unable to propagate logger
class FHIAimsMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class FHIAimsOutMappingParser(TextMappingParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    _gw_flag_map = {
        'gw': 'G0W0',
        'gw_expt': 'G0W0',
        'ev_scgw0': 'ev-scGW',
        'ev_scgw': 'ev-scGW',
        'scgw': 'scGW',
    }

    # FHI-aims XC control description -> standard functional name. The
    # `nomad-simulations` schema expands the name into LibXC components
    # (family/kind) and derives `jacobs_ladder`; the parser deliberately does
    # not resolve LibXC labels itself. Names normalize case-insensitively, so
    # the human-readable spelling used here resolves to the schema alias.
    _xc_name_map = {
        'Perdew-Wang parametrisation of Ceperley-Alder LDA': 'LDA',
        'Perdew-Zunger parametrisation of Ceperley-Alder LDA': 'PZ81',
        'VWN-LDA parametrisation of VWN5 form': 'VWN',
        'VWN-LDA parametrisation of VWN-RPA form': 'VWN-RPA',
        'AM05 gradient-corrected functionals': 'AM05',
        'BLYP functional': 'BLYP',
        'PBE gradient-corrected functionals': 'PBE',
        'PBEint gradient-corrected functional': 'PBEint',
        'PBEsol gradient-corrected functionals': 'PBEsol',
        'RPBE gradient-corrected functionals': 'RPBE',
        'revPBE gradient-corrected functionals': 'revPBE',
        'PW91 gradient-corrected functionals': 'PW91',
        'M06-L gradient-corrected functionals': 'M06-L',
        'M11-L gradient-corrected functionals': 'M11-L',
        'TPSS gradient-corrected functionals': 'TPSS',
        'TPSSloc gradient-corrected functionals': 'TPSSloc',
        'hybrid B3LYP functional': 'B3LYP',
        'HSE': 'HSE03',
        'HSE-functional': 'HSE06',
        'hybrid-PBE0 functionals': 'PBE0',
        'hybrid-PBEsol0 functionals': 'PBEsol0',
        'Hybrid M06 gradient-corrected functionals': 'M06',
        'Hybrid M06-2X gradient-corrected functionals': 'M06-2X',
        'Hybrid M06-HF gradient-corrected functionals': 'M06-HF',
        'Hybrid M08-HX gradient-corrected functionals': 'M08-HX',
        'Hybrid M08-SO gradient-corrected functionals': 'M08-SO',
        'Hybrid M11 gradient-corrected functionals': 'M11',
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

    def get_functional_key(self, xc: str) -> str | None:
        """Standard functional name for this FHI-aims XC control string. The
        schema expands the name into LibXC components (family/kind) and derives
        `jacobs_ladder`; returns `None` for an unmapped string."""
        return self._xc_name_map.get(xc)

    def get_periodic_boundary_conditions(self, source: dict[str, Any]) -> list[bool]:
        return [source.get('lattice_vectors') is not None] * 3

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

    def get_band_structures(
        self, source: list[dict[str, Any]], params: dict[str, Any]
    ) -> list[dict[str, Any]]:
        band_structures = []
        for spin_channel, eig in enumerate(self.get_eigenvalues(source, params)):
            values = eig.get('eigenvalues')
            if values is None:
                continue
            band_structures.append(
                dict(
                    value=values,
                    occupation=eig.get('occupations'),
                    spin_channel=spin_channel,
                )
            )
        return band_structures

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

    def get_k_offset_with_default(self, k_offset: np.ndarray | None) -> np.ndarray:
        """
        Return k_offset or FHI-aims default [0,0,0] (Gamma-centered).

        FHI-aims uses Gamma-centered grids by default when k_offset
        is not specified in control.in.
        """
        if k_offset is None:
            return np.array([0.0, 0.0, 0.0])
        return k_offset

    def get_all_criteria(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        criteria = []
        for fcriteria in [
            self.get_scf_energy_criterion,
            self.get_scf_density_criterion,
            self.get_scf_eigenvalues_criterion,
        ]:
            data = fcriteria(source)
            if data:
                criteria.append(data)
        if not criteria:
            max_iterations = source.get('max_scf_iterations')
            if max_iterations is not None:
                criteria.append({'n_max_iterations': max_iterations})
        return criteria

    def get_scf_energy_criterion(self, source: dict[str, Any]) -> dict[str, Any]:
        """
        Extract energy convergence criterion for SCF.

        Returns:
            Dictionary representing energy SelfConsistency section, or empty dict.
        """
        conv_energy = source.get('convergence_energy')
        if conv_energy is None:
            return {}

        result = {
            'name': 'total_energy_change',
            'threshold_change': conv_energy,
        }

        max_iterations = source.get('max_scf_iterations')
        if max_iterations is not None:
            result['n_max_iterations'] = max_iterations

        return result

    def get_scf_density_criterion(self, source: dict[str, Any]) -> dict[str, Any]:
        """
        Extract density convergence criterion for SCF.

        Returns:
            Dictionary representing density SelfConsistency section, or empty dict.
        """
        conv_density = source.get('convergence_density')
        if conv_density is None:
            return {}

        result = {
            'name': 'charge_density_change',
            'threshold_change': conv_density,
            'threshold_change_unit': 'dimensionless',
        }

        max_iterations = source.get('max_scf_iterations')
        if max_iterations is not None:
            result['n_max_iterations'] = max_iterations

        return result

    def get_scf_eigenvalues_criterion(self, source: dict[str, Any]) -> dict[str, Any]:
        """
        Extract eigenvalues convergence criterion for SCF.

        Returns:
            Dictionary representing eigenvalues SelfConsistency section, or empty dict.
        """
        conv_eigenvalues = source.get('convergence_eigenvalues')
        if conv_eigenvalues is None:
            return {}

        result = {
            'name': 'sum_eigenvalues_change',
            'threshold_change': conv_eigenvalues,
        }

        max_iterations = source.get('max_scf_iterations')
        if max_iterations is not None:
            result['n_max_iterations'] = max_iterations

        return result

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        scf_iterations = source.get('self_consistency', [])
        if not scf_iterations:
            return {}

        delta_energies_total = []
        delta_charge_abs = []
        durations = []
        delta_sum_eigenvalues = []

        for iteration in scf_iterations:
            convergence = iteration.get('scf_convergence', {})
            delta_energy = convergence.get('Change of total energy')
            if delta_energy is not None:
                delta_energies_total.append(abs(delta_energy))

            delta_density = convergence.get('Change of charge density')
            if delta_density is not None:
                delta_charge_abs.append(
                    abs(float(delta_density)) * ureg.elementary_charge
                )

            delta_eig = convergence.get('Change of sum of eigenvalues')
            if delta_eig is not None:
                delta_sum_eigenvalues.append(delta_eig)

            duration = iteration.get('time_calculation')
            if duration is not None:
                durations.append(float(duration))

        scf_steps = {}
        if delta_energies_total:
            scf_steps['delta_energies_total'] = delta_energies_total
        if delta_charge_abs:
            scf_steps['delta_charge_abs'] = delta_charge_abs
        if len(durations) == len(scf_iterations):
            scf_steps['durations'] = durations
        if delta_sum_eigenvalues:
            scf_steps['code_specific_quantities'] = {
                'delta_sum_eigenvalues': delta_sum_eigenvalues
            }
        return scf_steps

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
    annotation_key: str = fhiaims.TEXT_KEY
    geometry_parser = GeometryParser()
    control_parser = ControlParser()

    def build_phonopy_object(
        self, tolerance=1e-6
    ) -> tuple[Phonopy, dict[str, EntryArchive]]:
        """
        Generate Phonopy object from FHI-aims files, also returns the archives
        containing the force sets.
        """

        def is_equal(reference: Atoms, calculated: Atoms) -> bool:
            """
            Compare ase Atoms objects.
            """
            if len(reference) != len(calculated):
                return False
            if (
                reference.get_atomic_numbers() != calculated.get_atomic_numbers()
            ).any():
                return False
            if (abs(reference.get_cell() - calculated.get_cell()) > tolerance).any():
                return False
            # get normalized positions, wrapped to the bounding cell
            ref_pos = reference.get_scaled_positions() % 1.0
            cal_pos = calculated.get_scaled_positions() % 1.0
            # resolve coordinates at the boundary
            ref_pos = np.where(ref_pos != 1.0, ref_pos, 0.0)
            cal_pos = np.where(cal_pos != 1.0, cal_pos, 0.0)
            if (abs(ref_pos - cal_pos) > tolerance).any():
                return False
            return True

        def get_forces(
            mainfile: str, supercell: PhonopyAtoms, archive: EntryArchive = None
        ) -> tuple[EntryArchive, np.ndarray]:
            """
            Load archive in upload corresponding to mainfile.
            """
            if archive is None:
                archive = self.archive.m_context.resolve_archive(
                    f'../upload/archive/mainfile/{mainfile}'
                )
            # check if supercell match calculation cell
            calc_cell: Atoms = archive.data.model_system[-1].to_ase_atoms(
                logger=self.logger
            )
            supercell_atoms = Atoms(
                positions=supercell.positions,
                cell=supercell.cell,
                symbols=supercell.symbols,
                pbc=True,
            )

            if not is_equal(supercell_atoms, calc_cell):
                self.logger.error('Phonopy supercell does not match calculation.')

            forces = (
                archive.data.outputs[-1]
                .total_forces[-1]
                .value.to('eV/angstrom')
                .magnitude
            )

            return archive, forces

        maindir = os.path.dirname(os.path.dirname(self.mainfile))

        self.control_parser.mainfile = os.path.join(maindir, 'control.in')
        supercell_matrix = self.control_parser.get('supercell')
        displacement = self.control_parser.get('displacement', 0.001)
        sym = self.control_parser.get('symmetry_thresh', 1e-6)

        self.geometry_parser.mainfile = os.path.join(maindir, 'geometry.in')
        unit_atoms = PhonopyAtoms(atoms=self.geometry_parser.get_atoms())

        phonopy_obj = Phonopy(
            unit_atoms, supercell_matrix, symprec=sym, calculator='fhi-aims'
        )
        phonopy_obj.generate_displacements(distance=displacement)
        supercells = phonopy_obj.supercells_with_displacements

        force_sets = []
        n_pad = int(np.ceil(np.log10(len(supercells) + 1))) + 1
        force_archives: dict[str, EntryArchive] = {}
        supercell = supercells[0]
        for n in range(3):
            # for n, supercell in enumerate(supercells):
            calc_dir = f'phonopy-FHI-aims-displacement-{str(n + 1).zfill(n_pad)}'
            calc_file = search_files(
                '*out', os.path.join(maindir, calc_dir), re_pattern=f'{calc_dir}.out'
            )
            if not calc_file:
                self.logger.error('No FHI-aims phonon calculation file found.')
                break

            match = re.match(r'.+?/raw/(.+)', calc_file[0])
            if not match:
                break

            # TODO put a try here, wait until archive is processed
            archive, forces = get_forces(
                match.group(1), supercell, self.archive if n == 0 else None
            )
            force_archives[f'supercell {n}'] = archive
            force_sets.append(forces)

        try:
            phonopy_obj.forces = force_sets
            phonopy_obj.produce_force_constants()
        except Exception:
            self.logger.error('Error producing force constants.')

        return phonopy_obj, force_archives

    def write_to_archive(  # noqa: PLR0915, PLR0912
        self,
    ) -> None:
        out_parser = FHIAimsOutMappingParser()
        out_parser.text_parser = FHIAimsOutFileParser()
        out_parser.filepath = self.mainfile

        archive_handler = FHIAimsMetainfoParser()
        archive_handler.annotation_key = self.annotation_key

        self.archive.data = Simulation(program=Program(name='FHI-aims'))
        archive_handler.data_object = self.archive.data

        out_parser.convert(archive_handler, remove=False)

        # separate parsing of dos due to a problem with mapping physical
        # property variables
        archive_handler.annotation_key = fhiaims.TEXT_DOS_KEY
        out_parser.convert(archive_handler, remove=False)

        # workflow
        if out_parser.data.get('geometry_optimization'):
            workflow_key = 'geo_opt_workflow'
            self.archive.workflow2 = GeometryOptimization()
            self.archive.workflow2.method = GeometryOptimizationMethod()
            self.archive.workflow2.method.optimization_method = out_parser.data.get(
                'geometry_relaxation_method'
            )
            force_threshold = out_parser.data.get('convergence_forces')
            if force_threshold is not None:
                self.archive.workflow2.method.convergence_targets = [
                    ForceConvergenceTarget(
                        threshold=force_threshold,
                        threshold_type='maximum',
                    )
                ]
            energy_threshold = out_parser.data.get('convergence_energy')
            if energy_threshold is not None:
                # Handle both pint Quantity (with units) and plain float
                if hasattr(energy_threshold, 'units'):
                    threshold_value = energy_threshold
                else:
                    threshold_value = energy_threshold * ureg.eV

                self.archive.workflow2.method.single_point_convergence_targets = [
                    EnergyConvergenceTarget(
                        threshold=threshold_value,
                        threshold_type='absolute',
                    )
                ]
        elif out_parser.data.get('molecular_dynamics'):
            workflow_key = 'md_workflow'
            self.archive.workflow2 = MolecularDynamics()
        else:
            workflow_key = None
            self.archive.workflow2 = SinglePoint()
            self.archive.workflow2.method = SinglePointMethod()
            energy_threshold = out_parser.data.get('convergence_energy')
            if energy_threshold is not None:
                self.archive.workflow2.method.convergence_targets = [
                    EnergyConvergenceTarget(
                        threshold=energy_threshold,
                        threshold_type='absolute',
                    )
                ]
        if workflow_key:
            archive_handler.data_object = self.archive.workflow2
            archive_handler.annotation_key = workflow_key
            out_parser.convert(archive_handler)

        gw_archive = self.child_archives.get('GW') if self.child_archives else None
        if gw_archive is not None:
            self.archive.workflow2.name = 'DFT'

            # GW single point
            parser = FHIAimsArchiveWriter()
            parser.annotation_key = fhiaims.TEXT_GW_KEY
            parser.write(self.mainfile, gw_archive, self.logger)
            gw_archive.workflow2.name = 'GW'

            # DFT-GW workflow
            gw_workflow_archive = self.child_archives.get('GW_workflow')
            gw_workflow_archive.workflow2 = DFTGWWorkflow(
                tasks=[
                    TaskReference(task=self.archive.workflow2),
                    TaskReference(task=gw_archive.workflow2),
                ]
            )

        phonon_archive = (
            self.child_archives.get('phonon') if self.child_archives else None
        )
        if phonon_archive is not None:
            # TODO generalize this to a parser/normalizer
            # create phonony object
            phonopy_obj, force_archives = self.build_phonopy_object()
            # run phonopy calculation and fill phonon archive
            try:
                phonopy_obj_to_archive(phonopy_obj, phonon_archive, self.logger)
            except Exception:
                self.logger.error('Failed to run phonopy.')
            # link calculations and phonon workflow
            phonon_workflow_archive = self.child_archives.get('phonon_workflow')
            phonon_workflow_archive.workflow2 = Phonon(
                # unit cell
                inputs=[
                    Link(
                        name='Input system', section=phonon_archive.data.model_system[0]
                    )
                ],
                tasks=[
                    # SimulationTask(name='Supercell generation'),
                    *[
                        SimulationTaskReference(
                            task=a.workflow2,
                        )
                        for k, a in force_archives.items()
                        if a.workflow2
                    ],
                    SimulationTaskReference(task=phonon_archive.workflow2),
                ],
            )
            phonon_workflow_archive.workflow2.tasks[0].inputs.extend(
                phonon_workflow_archive.workflow2.inputs
            )
            for archive in force_archives.values():
                if not archive.workflow2:
                    continue
                archive.workflow2.normalize(archive, self.logger)
                phonon_workflow_archive.workflow2.tasks[0].outputs.extend(
                    archive.workflow2.inputs
                )

        # close file contexts
        out_parser.close()
        archive_handler.close()


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
    ) -> bool | Iterable[str]:
        is_mainfile = super().is_mainfile(
            filename=filename,
            mime=mime,
            buffer=buffer,
            decoded_buffer=decoded_buffer,
            compression=compression,
        )
        children = []
        if is_mainfile:
            # gw calculation
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
                children = ['GW', 'GW_workflow']
            if not children:
                # phonon calculation
                match = re.search(r'.*/phonopy-FHI-aims-displacement-0+1/.*', filename)
                if match:
                    children = ['phonon', 'phonon_workflow']

        if children:
            # TODO not possible at the moment to redefine level
            # self.level = 1
            self.creates_children = True

        return children or is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, EntryArchive] = None,
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
