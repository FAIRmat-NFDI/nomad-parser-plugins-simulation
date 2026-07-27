import os
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.workflow import (
    GeometryOptimization,
    MolecularDynamics,
    SinglePoint,
)
from nomad_simulations.schema_packages.workflow.general import (
    EnergyConvergenceTarget,
    ForceConvergenceTarget,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.single_point import SinglePointMethod

from nomad_simulation_parsers.parsers.utils.general import (
    as_list,
    calculate_band_gap_from_occupations,
)
from nomad_simulation_parsers.schema_packages import vasp

from .common import get_functional_key as _functional_key_from_params
from .outcar_parser import OutcarArchiveWriter

LOGGER = get_logger(__name__)
N_SPIN_CHANNELS = 2
EIGENVALUE_COMPONENTS = 2


# TODO temporary fix for structlog unable to propagate logger
class VASPMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class VasprunParser(XMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def mix_alpha(self, mix: float, cond: bool) -> float:
        return mix if cond else 0

    def _extract_spin_sets(
        self, source: dict[str, Any] | None
    ) -> list[list[dict[str, Any]]]:
        if source is None:
            return []
        node = source.get('array', source)
        if isinstance(node, dict):
            node = node.get('set')
        top_level = as_list(node)
        if not top_level:
            return []

        # XML payload can contain one or more wrapper `set` containers.
        while (
            len(top_level) == 1
            and isinstance(top_level[0], dict)
            and top_level[0].get('r') is None
            and top_level[0].get('set') is not None
        ):
            top_level = as_list(top_level[0].get('set'))
            if not top_level:
                return []

        # Non-spin structure: top-level is already a list of k-point entries with `r`.
        if all(
            isinstance(item, dict) and item.get('r') is not None for item in top_level
        ):
            return [top_level]

        # Spin structure: top-level contains one container per spin channel.
        spin_sets = []
        for item in top_level:
            if not isinstance(item, dict):
                continue
            kpoint_sets = as_list(item.get('set'))
            if kpoint_sets and all(
                isinstance(kpt, dict) and kpt.get('r') is not None
                for kpt in kpoint_sets
            ):
                spin_sets.append(kpoint_sets)
        return spin_sets

    def get_eigenvalues(self, source: dict[str, Any] | None) -> list[dict[str, Any]]:
        spin_sets = self._extract_spin_sets(source)
        if not spin_sets:
            return []

        eigenvalue_array_ndim = 3  # (n_kpoints, n_bands, components)
        eigenvalues = []
        is_spin_polarized = len(spin_sets) == N_SPIN_CHANNELS
        for spin_channel, kpoint_sets in enumerate(spin_sets):
            rows = [kpt.get('r') for kpt in kpoint_sets if isinstance(kpt, dict)]
            if not rows:
                continue
            try:
                data = np.asarray(rows, dtype=float)
            except Exception:
                continue
            if (
                data.ndim != eigenvalue_array_ndim
                or data.shape[2] < EIGENVALUE_COMPONENTS
            ):
                continue

            eigs = data[:, :, 0]
            occs = data[:, :, 1]
            entry = dict(
                value=eigs,
                occupation=occs,
                n_levels=eigs.shape[1],
            )
            if is_spin_polarized:
                entry['spin_channel'] = spin_channel
            eigenvalues.append(entry)
        return eigenvalues

    def get_band_gaps(self, source: dict[str, Any] | None) -> list[dict[str, Any]]:
        """Calculate band gaps from eigenvalues using common utility."""
        result = []
        for eigenvalues in self.get_eigenvalues(source):
            eigs = eigenvalues.get('value')
            occs = eigenvalues.get('occupation')
            spin_channel = eigenvalues.get('spin_channel')

            # Use common utility for band gap calculation
            gap_result = calculate_band_gap_from_occupations(
                eigs, occs, spin_channel=spin_channel
            )
            if gap_result is not None:
                result.append(gap_result)

        return result

    def get_total_dos(self, source: dict[str, Any] | None) -> list[dict[str, Any]]:  # noqa: PLR0912
        def get_fermi_energy(source: dict[str, Any] | None) -> float | None:
            if not isinstance(source, dict):
                return None
            fermi = source.get('i')
            if isinstance(fermi, dict):
                value = fermi.get(self.value_key)
                try:
                    return float(value) if value is not None else None
                except Exception:
                    return None
            for item in as_list(fermi):
                if not isinstance(item, dict):
                    continue
                if item.get(f'{self.attribute_prefix}name') != 'efermi':
                    continue
                value = item.get(self.value_key)
                try:
                    return float(value) if value is not None else None
                except Exception:
                    return None
            return None

        dos = source
        if not isinstance(dos, dict):
            return []
        total = dos.get('total')
        if not isinstance(total, dict):
            return []

        dos_array_ndim = 2  # (n_points, columns)
        node = total.get('array', total)
        if isinstance(node, dict):
            node = node.get('set')
        spin_entries = as_list(node)
        while (
            len(spin_entries) == 1
            and isinstance(spin_entries[0], dict)
            and spin_entries[0].get('r') is None
            and spin_entries[0].get('set') is not None
        ):
            spin_entries = as_list(spin_entries[0].get('set'))
            if not spin_entries:
                return []
        if not all(
            isinstance(entry, dict) and entry.get('r') is not None
            for entry in spin_entries
        ):
            return []

        efermi = get_fermi_energy(dos)
        is_spin_polarized = len(spin_entries) == N_SPIN_CHANNELS
        dos_sections = []
        for spin_channel, dos_entry in enumerate(spin_entries):
            rows = dos_entry.get('r')
            if rows is None:
                continue
            try:
                data = np.asarray(rows, dtype=float)
            except Exception:
                continue
            if data.ndim != dos_array_ndim or data.shape[1] < EIGENVALUE_COMPONENTS:
                continue
            # VASP may print spin-down DOS with negative sign for plotting.
            # The schema expects non-negative DOS intensities.
            entry = dict(energies=data[:, 0], value=np.abs(data[:, 1]))
            if is_spin_polarized:
                entry['spin_channel'] = spin_channel
            if efermi is not None:
                entry['energy_fermi'] = efermi
            dos_sections.append(entry)
        return dos_sections

    def get_energy_contributions(
        self, source: list[dict[str, Any]], **kwargs
    ) -> list[dict[str, Any]]:
        return [
            c
            for c in source
            if c.get(f'{self.attribute_prefix}name') not in kwargs.get('exclude', [])
        ]

    def get_data(self, source: dict[str, Any], **kwargs) -> Any:
        if source.get(self.value_key):
            return source[self.value_key]
        path = kwargs.get('path')
        if path is None:
            return

        parser = Path(path=path)
        return parser.get_data(source)

    def get_forces(self, source: dict[str, Any]) -> dict[str, Any]:
        value = self.get_data(source, path='.varray.v')
        if value is None:
            return {}
        return dict(forces=value, npoints=len(value), rank=[3])

    def reshape_array(self, source: np.ndarray, shape_rest: tuple = (3,)) -> np.ndarray:
        if source is None:
            return
        return np.reshape(
            source, (np.size(source) // int(np.prod(shape_rest)), *shape_rest)
        )

    def _get_parameter(self, name: str, section: str | None = None) -> Any:
        parameters = self.data.get('modeling', {}).get('parameters', {})
        separators = parameters.get('separator', [])
        if not isinstance(separators, list):
            return None
        for separator in separators:
            if section and separator.get(self.attribute_prefix + 'name') != section:
                continue
            for quantity in separator.get('i', []):
                if quantity.get(self.attribute_prefix + 'name') == name:
                    return quantity.get(self.value_key)
        return None

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:  # noqa: PLR0912
        scf_iterations = source.get('scstep', [])
        if not scf_iterations:
            return {}

        energies_total = []
        durations = []
        for step in scf_iterations:
            entries = step.get('energy', {}).get('i', [])
            energy_total = None
            for entry in entries:
                if entry.get(self.attribute_prefix + 'name') == 'e_fr_energy':
                    energy_total = entry.get(self.value_key)
                    break
            if energy_total is not None:
                energies_total.append(float(energy_total) * ureg.eV)

            total_time = None
            for timing in step.get('time', []):
                if timing.get(self.attribute_prefix + 'name') != 'total':
                    continue
                value = timing.get(self.value_key)
                if isinstance(value, np.ndarray | list | tuple):
                    if len(value) > 1:
                        total_time = float(value[1])
                    elif len(value) == 1:
                        total_time = float(value[0])
                elif isinstance(value, int | float):
                    total_time = float(value)
                break
            if total_time is not None:
                durations.append(total_time)

        if not energies_total:
            return {}

        scf_steps = {'energies_total': energies_total}
        if len(energies_total) > 1:
            delta_energies_total = []
            for idx in range(1, len(energies_total)):
                delta_energies_total.append(
                    abs(energies_total[idx] - energies_total[idx - 1])
                )
            scf_steps['delta_energies_total'] = delta_energies_total
        if len(durations) == len(energies_total):
            scf_steps['durations'] = durations
        return scf_steps

    def build_workflow(self):
        nsw = self._get_parameter('NSW', section='ionic')
        ibrion = -1 if nsw == 0 else self._get_parameter('IBRION', section='ionic')
        if ibrion is None:
            ibrion = -1

        ediff = self._get_parameter('EDIFF', section='electronic')

        def single_point_convergence() -> list[EnergyConvergenceTarget]:
            if ediff is None:
                return []
            return [
                EnergyConvergenceTarget(
                    threshold=float(ediff) * ureg.eV,
                    threshold_type='absolute',
                )
            ]

        if int(ibrion) == -1:
            workflow = SinglePoint()
            workflow.method = SinglePointMethod()
            convergence = single_point_convergence()
            if convergence:
                workflow.method.convergence_targets = convergence
            return workflow
        if int(ibrion) == 0:
            return MolecularDynamics()

        workflow = GeometryOptimization()
        workflow.method = GeometryOptimizationMethod()

        ediffg = self._get_parameter('EDIFFG', section='ionic')
        convergence_targets = []
        if ediffg is not None:
            ediffg = float(ediffg)
            if ediffg > 0:
                convergence_targets.append(
                    EnergyConvergenceTarget(
                        threshold=ediffg * ureg.eV,
                        threshold_type='absolute',
                    )
                )
            else:
                convergence_targets.append(
                    ForceConvergenceTarget(
                        threshold=abs(ediffg) * ureg.eV / ureg.angstrom,
                        threshold_type='maximum',
                    )
                )
        if convergence_targets:
            workflow.method.convergence_targets = convergence_targets

        sp_convergence = single_point_convergence()
        if sp_convergence:
            workflow.method.single_point_convergence_targets = sp_convergence
        return workflow

    def get_atoms(self, arrays: Any = None) -> list[dict[str, str]]:
        arrays = as_list(arrays)
        atoms_array = next(
            (
                array
                for array in arrays
                if array.get(f'{self.attribute_prefix}name') == 'atoms'
            ),
            None,
        )
        if atoms_array is None:
            return []
        rows = atoms_array.get('set', {}).get('rc', [])
        if isinstance(rows, dict):
            rows = [rows]
        atoms = []
        for row in rows:
            values = row.get('c')
            if isinstance(values, list) and len(values) > 0:
                symbol = values[0]
            else:
                symbol = values
            if symbol is None:
                continue
            atoms.append({'label': str(symbol).strip()})
        return atoms

    def get_positions(
        self, positions: Any = None, lattice_vectors: Any = None
    ) -> np.ndarray | None:
        if positions is None:
            return None
        if lattice_vectors is None:
            return positions
        return np.dot(np.asarray(positions), np.asarray(lattice_vectors))

    def get_periodic_boundary_conditions(self) -> list[bool]:
        return [True, True, True]

    def _find_parameter(self, name: str) -> Any:
        """Recursively search the vasprun `<parameters>` tree for an `<i name=...>`,
        returning its value typed by the element's `type` attribute (`logical` ->
        bool, `int`/`float` -> number), mirroring the OUTCAR parameter typing.

        The XC-relevant tags (`GGA`, `METAGGA`, `LHFCALC`, `AEXX`, ...) live in
        different, sometimes nested, `<separator>` blocks, so a flat lookup of a
        single bound node is not enough.
        """

        def typed(item: dict[str, Any]) -> Any:
            value = item.get(self.value_key)
            if not isinstance(value, str):
                return value
            text = value.strip()
            match item.get(f'{self.attribute_prefix}type'):
                case 'logical':
                    return text.strip('.').upper() in ('T', 'TRUE')
                case 'string':
                    return text
                case 'int':
                    return int(text) if text.lstrip('-').isdigit() else None
                case _:
                    # untyped numeric `<i>` default to float in vasprun
                    try:
                        return float(text)
                    except ValueError:
                        return text

        def search(node: Any) -> Any:
            if not isinstance(node, dict):
                return None
            for item in as_list(node.get('i')):
                if (
                    isinstance(item, dict)
                    and item.get(f'{self.attribute_prefix}name') == name
                ):
                    return typed(item)
            for sub in as_list(node.get('separator')):
                found = search(sub)
                if found is not None:
                    return found
            return None

        parameters = self.data.get('modeling', {}).get('parameters', {})
        return search(parameters)

    def get_functional_key(self, source: Any = None) -> str | None:
        params = {
            key: self._find_parameter(key)
            for key in (
                'GGA',
                'METAGGA',
                'LHFCALC',
                'AEXX',
                'AGGAC',
                'ALDAC',
                'HFSCREEN',
            )
        }
        return _functional_key_from_params(params)


class XMLArchiveWriter(ArchiveWriter):
    def _has_electronic_outputs(self) -> bool:
        outputs = getattr(getattr(self.archive, 'data', None), 'outputs', None) or []
        for output in outputs:
            if getattr(output, 'electronic_band_structures', None):
                return True
            if getattr(output, 'electronic_dos', None):
                return True
            if getattr(output, 'electronic_band_gaps', None):
                return True
        return False

    # TODO(mapping-migration): replace this XML->OUTCAR backfill with a
    # mapping-driven source merge when XML fixtures with missing electronic
    # payloads are supported in mappings; see the tracking issue. Disabling
    # the backfill regresses tests/parsers/test_vasp_parser.py::
    # test_vasprun_backfills_electronic_outputs_from_outcar_when_xml_missing.
    def _backfill_from_outcar(self) -> None:
        if self._has_electronic_outputs():
            return

        outcar_path = os.path.join(os.path.dirname(self.mainfile), 'OUTCAR')
        if not os.path.isfile(outcar_path):
            return

        outcar_archive = EntryArchive()
        OutcarArchiveWriter().write(
            outcar_path, outcar_archive, self.logger, self.child_archives
        )

        outcar_outputs = (
            getattr(getattr(outcar_archive, 'data', None), 'outputs', None) or []
        )
        if not outcar_outputs:
            return

        outcar_output = outcar_outputs[0]
        target_outputs = (
            getattr(getattr(self.archive, 'data', None), 'outputs', None) or []
        )
        if not target_outputs:
            self.archive.data.outputs = [outcar_output]
            return

        target_output = target_outputs[0]
        for quantity_name in (
            'electronic_band_structures',
            'electronic_band_gaps',
            'electronic_dos',
        ):
            if getattr(target_output, quantity_name, None):
                continue
            value = getattr(outcar_output, quantity_name, None)
            if value:
                setattr(target_output, quantity_name, value)

    def write_to_archive(self) -> None:
        data_parser = VASPMetainfoParser()
        data_parser.data_object = Simulation()

        xml_parser = VasprunParser(filepath=self.mainfile)

        data_parser.annotation_key = vasp.XML_KEY
        xml_parser.convert(data_parser)

        data_parser.annotation_key = vasp.XML2_KEY
        xml_parser.convert(data_parser)

        self.archive.data = data_parser.data_object
        self.archive.workflow2 = xml_parser.build_workflow()

        self._backfill_from_outcar()

        # close file objects
        data_parser.close()
        xml_parser.close()
