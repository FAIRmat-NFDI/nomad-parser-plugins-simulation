from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    pass

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

from nomad_simulation_parsers.schema_packages import vasp

LOGGER = get_logger(__name__)


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

    def get_eigenvalues(self, array: list) -> dict[str, Any]:
        if array is None:
            return {}
        transposed = np.transpose(array)
        return dict(eigenvalues=transposed[0], occupations=transposed[1])

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
