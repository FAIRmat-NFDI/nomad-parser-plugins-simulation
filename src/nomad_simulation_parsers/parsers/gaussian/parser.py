from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from importlib import reload

import numpy as np
from ase.data import chemical_symbols
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.schema_packages import gaussian
from nomad_simulation_parsers.schema_packages.utils import remove_mapping_annotations

from .text_parser import GaussianOutReader

LOGGER = get_logger(__name__)


class OutParser(MappingTextParser):
    """
    Couples GaussianOutReader (regex) with a few convenience getters that
    the mapping rules will call.
    """

    def __init__(self):
        super().__init__(text_parser=GaussianOutReader())

    def get_program_data(self, src: dict[str, Any]) -> dict[str, Any]:
        program = src.get('program')
        if isinstance(program, list):
            if program and isinstance(program[0], list):
                program = program[0]
            version = ' '.join([str(val) for val in program if val])
        elif program is not None:
            version = str(program)
        else:
            version = None
        return {
            'name': 'Gaussian',
            'version': version,
        }

    @staticmethod
    def _get_last(values):
        if values is None:
            return None
        if isinstance(values, list):
            return values[-1] if values else None
        return values

    @staticmethod
    def _get_total_spin(multiplicity: Any) -> Any:
        if multiplicity is None:
            return None
        try:
            return int(multiplicity) - 1
        except Exception:
            try:
                return float(multiplicity) - 1
            except Exception:
                return None

    @staticmethod
    def _select_orientation(system: dict[str, Any]) -> Any:
        return (
            system.get('standard_orientation')
            or system.get('input_orientation')
            or system.get('z_matrix_orientation')
        )

    def get_atoms(self, src):
        runs = src.get('run') or []
        runs = runs if isinstance(runs, list) else [runs]
        systems = []

        for run in runs:
            for system in run.get('system', []) or []:
                orientation = self._select_orientation(system)
                if orientation is None:
                    continue
                arr = np.asarray(orientation, dtype=float)
                if arr.size == 0:
                    continue

                positions = arr[:, -3:]
                atomic_numbers = arr[:, 1].astype(int) if arr.shape[1] > 1 else []
                atoms = [{'atomic_number': int(z)} for z in atomic_numbers]
                systems.append({'positions': positions, 'particle_states': atoms})

        return systems

    def get_outputs(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        runs = src.get('run') or []
        runs = runs if isinstance(runs, list) else [runs]
        outputs: list[dict[str, Any]] = []

        for run in runs:
            for system in run.get('system', []) or []:
                output: dict[str, Any] = {}

                energy = self._get_last(system.get('energy_total'))
                if energy is not None:
                    output['total_energies'] = [dict(value=energy)]

                forces = system.get('forces')
                if forces is not None:
                    output['total_forces'] = [dict(value=forces)]

                if output:
                    outputs.append(output)

        return outputs


class GaussianParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] | None = None,
    ) -> None:
        
        # Clean up any mapping annotations that might have been left by other parsers
        remove_mapping_annotations(Simulation.m_def)
        reload(gaussian)

        reader = OutParser()
        reader.filepath = mainfile

        meta = MetainfoParser(data_object=Simulation())
        meta.annotation_key = 'out'
        # Need to traverse Simulation -> ModelSystem -> particle_states
        # Depth: Simulation (0) -> model_system (1) -> particle_states (2) -> AtomsState quantities (3)
        meta.max_nested_level = 3

        reader.convert(meta)
        archive.data = meta.data_object

        meta.close()
        reader.close()

 
        # Remove ORCA mapping annotations to avoid interfering with other parsers.
        remove_mapping_annotations(Simulation.m_def)       
