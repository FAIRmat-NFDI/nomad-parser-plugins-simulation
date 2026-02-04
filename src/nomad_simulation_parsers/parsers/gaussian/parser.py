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
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_system import ModelSystem
from nomad_simulations.schema_packages import model_method
from nomad_simulations.schema_packages.atoms_state import AtomsState
from nomad.units import ureg

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
        # Parse everything; helper getters rely on full data presence.
        self.parse_only_required = False
        self.text_parser.parse_only_required = False
        self.text_parser.findlazy = False

    def load_file(self):
        if self.filepath:
            self.text_parser.findlazy = False
            self.text_parser.mainfile = self.filepath
        return self.text_parser

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

    def get_dft_method(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        """
        Parse Gaussian route section and return a list of DFT ModelMethod objects.
        Skips non-DFT methods (CC, MPn, Gx, HF, etc.).
        """
        runs = src.get('run') or []
        runs = runs if isinstance(runs, list) else [runs]
        methods: list[dict[str, Any]] = []

        skip = {
            'N',
            'P',
            'OPT',
            'FREQ',
            'TD',
            'INT',
            'SCF',
            'POP',
            'CHARGE',
            'GEOM',
            'GUESS',
            'SCRF',
            'IOP',
            'UNITS',
            'GFINPUT',
            'GFPRINT',
            'SYMMETRY',
            'NOSYMM',
            'NOPRINT',
        }
        exclude_prefixes = (
            'CC',
            'QC',
            'MP',
            'G1',
            'G2',
            'G3',
            'G4',
            'CBS',
            'HF',
            'RHF',
            'UHF',
            'ROHF',
            'CI',
            'CIS',
        )
        dft_whitelist = {
            'LDA',
            'SVWN',
            'BLYP',
            'BP86',
            'B3LYP',
            'PBE',
            'PBEPBE',
            'PBE0',
            'PBE1PBE',
            'PW91',
            'HSEH1PBE',
            'HSE06',
            'LC-WPBE',
            'CAM-B3LYP',
            'M06',
            'M06L',
            'M062X',
            'M11',
            'M11L',
            'MN12L',
            'MN12SX',
            'N12',
            'N12SX',
            'TPSSTPSS',
            'TPSSH',
            'SCAN',
            'WB97X',
            'WB97XD',
            'B97D',
            'B97D3',
            'APFD',
        }

        for run in runs:
            route = run.get('x_gaussian_settings_corrected', '') or ''
            route_up = route.upper()

            functional = None
            for token in route_up.replace('#', ' ').split():
                if '=' in token:
                    continue
                core = token.split('/')[0]
                core = core.lstrip('RUO')
                if len(core) < 3 or core in skip:
                    continue
                if core.startswith(exclude_prefixes):
                    continue
                if core not in dft_whitelist:
                    continue
                functional = core
                break

            if functional is None:
                continue

            dft = model_method.DFT()
            dft.name = 'DFT'
            dft.xc = model_method.XCFunctional(functional_key=functional)
            methods.append(dft)

        return methods

    def _extract_energy(self, obj: Any) -> Any:
        if obj is None:
            return None
        get = getattr(obj, 'get', lambda k, d=None: d)
        ene = self._get_last(get('energy_total'))
        if ene is not None:
            return ene
        ene = self._get_last(get('energy_total_scf_iteration'))
        if ene is not None:
            return ene
        scf_final = get('scf_iteration_final')
        if scf_final is not None:
            ene = self._get_last(getattr(scf_final, 'get', lambda k, d=None: d)('energy_total'))
            if ene is not None:
                return ene
        return None

    def get_outputs(self, src: dict[str, Any]) -> list[dict[str, Any]]:
        runs = src.get('run') or []
        runs = runs if isinstance(runs, list) else [runs]
        outputs: list[dict[str, Any]] = []
        idx = 0
        for run in runs:
            for system in run.get('system', []) or []:
                output: dict[str, Any] = {}
                ene = self._extract_energy(system)
                if ene is not None:
                    energy_val = ene if hasattr(ene, 'magnitude') else ene * ureg.hartree
                    output['total_energies'] = [dict(value=energy_val)]
                forces = system.get('forces')
                if forces is None:
                    for calc in system.get('calculation', []) or []:
                        forces = calc.get('forces')
                        if forces is not None:
                            break
                if forces is not None:
                    output['total_forces'] = [dict(value=np.asarray(forces))]
                if output:
                    output['wall_start'] = float(idx)
                    output['wall_end'] = float(idx)
                    outputs.append(output)
                idx += 1
        return outputs

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
        for key in ('standard_orientation', 'input_orientation', 'z_matrix_orientation'):
            ori = system.get(key)
            if ori is None:
                continue
            # Avoid numpy truth-value ambiguity
            try:
                if hasattr(ori, 'size'):
                    if ori.size == 0:
                        continue
                elif hasattr(ori, '__len__') and len(ori) == 0:
                    continue
            except Exception:
                pass
            return ori
        return None

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
                atoms = []
                for z in atomic_numbers:
                    z_int = int(z)
                    atoms.append(AtomsState(atomic_number=z_int))
                systems.append({'positions': positions, 'particle_states': atoms})

        return systems


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
        # Traverse nested sections (model_system, outputs, particle_states, etc.)
        meta.max_nested_level = 3

        reader.convert(meta)
        archive.data = meta.data_object

        meta.close()
        reader.close()


        remove_mapping_annotations(Simulation.m_def)       
