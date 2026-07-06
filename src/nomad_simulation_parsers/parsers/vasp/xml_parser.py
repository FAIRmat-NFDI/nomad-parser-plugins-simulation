from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path, XMLParser
from nomad.units import ureg
from nomad.utils import get_logger
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
    calculate_band_gap_from_occupations,
)
from nomad_simulation_parsers.schema_packages import vasp

LOGGER = get_logger(__name__)
N_SPIN_CHANNELS = 2
EIGENVALUE_COMPONENTS = 2
EIGENVALUE_ARRAY_NDIM = 3
DOS_ARRAY_NDIM = 2


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

    def _as_list(self, value: Any) -> list[Any]:
        if value is None:
            return []
        if isinstance(value, list):
            return value
        return [value]

    def _extract_spin_sets(
        self, source: dict[str, Any] | None
    ) -> list[list[dict[str, Any]]]:
        if source is None:
            return []
        node = source.get('array', source)
        if isinstance(node, dict):
            node = node.get('set')
        top_level = self._as_list(node)
        if not top_level:
            return []

        # XML payload can contain one or more wrapper `set` containers.
        while (
            len(top_level) == 1
            and isinstance(top_level[0], dict)
            and top_level[0].get('r') is None
            and top_level[0].get('set') is not None
        ):
            top_level = self._as_list(top_level[0].get('set'))
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
            kpoint_sets = self._as_list(item.get('set'))
            if kpoint_sets and all(
                isinstance(kpt, dict) and kpt.get('r') is not None
                for kpt in kpoint_sets
            ):
                spin_sets.append(kpoint_sets)
        return spin_sets

    def _resolve_electronic_source(
        self, source: dict[str, Any] | None, key: str
    ) -> dict[str, Any] | None:
        if isinstance(source, dict) and source.get(key) is not None:
            return source
        calculations = self.data.get('modeling', {}).get('calculation', [])
        if not isinstance(calculations, list):
            return None
        for calculation in reversed(calculations):
            if isinstance(calculation, dict) and calculation.get(key) is not None:
                return calculation
        return None

    def get_eigenvalues(self, source: dict[str, Any] | None) -> list[dict[str, Any]]:
        source = self._resolve_electronic_source(source, 'eigenvalues')
        if source is None:
            return []
        source = source.get('eigenvalues')
        spin_sets = self._extract_spin_sets(source)
        if not spin_sets:
            return []

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
                data.ndim != EIGENVALUE_ARRAY_NDIM
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

    def _get_fermi_energy(self, source: dict[str, Any] | None) -> float | None:
        if not isinstance(source, dict):
            return None
        fermi = source.get('i')
        if isinstance(fermi, dict):
            value = fermi.get(self.value_key)
            try:
                return float(value) if value is not None else None
            except Exception:
                return None
        for item in self._as_list(fermi):
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

    def get_total_dos(self, source: dict[str, Any] | None) -> list[dict[str, Any]]:  # noqa: PLR0912
        source = self._resolve_electronic_source(source, 'dos')
        if not isinstance(source, dict):
            return []
        dos = source.get('dos')
        if not isinstance(dos, dict):
            return []
        total = dos.get('total')
        if not isinstance(total, dict):
            return []

        node = total.get('array', total)
        if isinstance(node, dict):
            node = node.get('set')
        spin_entries = self._as_list(node)
        while (
            len(spin_entries) == 1
            and isinstance(spin_entries[0], dict)
            and spin_entries[0].get('r') is None
            and spin_entries[0].get('set') is not None
        ):
            spin_entries = self._as_list(spin_entries[0].get('set'))
            if not spin_entries:
                return []
        if not all(
            isinstance(entry, dict) and entry.get('r') is not None
            for entry in spin_entries
        ):
            return []

        efermi = self._get_fermi_energy(dos)
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
            if data.ndim != DOS_ARRAY_NDIM or data.shape[1] < EIGENVALUE_COMPONENTS:
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

    def get_atoms(self) -> list[dict[str, str]]:
        modeling = self.data.get('modeling', {})
        atominfo = modeling.get('atominfo', {})
        arrays = atominfo.get('array', [])
        if isinstance(arrays, dict):
            arrays = [arrays]
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

    def get_positions(self, structure: dict[str, Any]) -> np.ndarray | None:
        positions = Path(path='.varray.v').get_data(structure)
        lattice_vectors = self.get_lattice_vectors(structure)
        if positions is None:
            return None
        if lattice_vectors is None:
            return positions
        return np.dot(np.asarray(positions), np.asarray(lattice_vectors))

    def get_lattice_vectors(self, structure: dict[str, Any]) -> np.ndarray | None:
        return Path(path='.crystal.varray[?"@name"==\'basis\'] | [0].v').get_data(
            structure
        )

    def get_periodic_boundary_conditions(self) -> list[bool]:
        return [True, True, True]

    def get_xc_functionals(self, source: dict[str, Any] | None) -> list[dict[str, str]]:
        if source is None:
            return []
        raw_params = source.get('i') if isinstance(source, dict) else None
        params = {}
        for item in raw_params or []:
            key = item.get(f'{self.attribute_prefix}name')
            value = item.get(self.value_key)
            if key is not None:
                params[key] = value
        if not params:
            return []
        # Reuse VASP OUTCAR XC mapping table to keep XML/OUTCAR behavior aligned.
        outcar_module = __import__(
            'nomad_simulation_parsers.parsers.vasp.outcar_parser',
            fromlist=['OutcarParser'],
        )
        return outcar_module.OutcarParser().get_xc_functionals(params)

    def get_ediff_unit(self) -> str:
        # VASP EDIFF is an energy threshold in eV.
        return 'electron_volt'


class XMLArchiveWriter(ArchiveWriter):
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

        # close file objects
        data_parser.close()
        xml_parser.close()
