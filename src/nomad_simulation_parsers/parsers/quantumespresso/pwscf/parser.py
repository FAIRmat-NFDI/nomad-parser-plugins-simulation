from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import pwscf

from ..parser import MainfileTextParser, MainfileXMLParser
from .file_parser import PWSCFFileParser

LOGGER = get_logger(__name__)


class PWSCFMainfileTextParser(MainfileTextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_force_contributions(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        keys = ['dispersion']
        return [
            dict(name=key, value=source[f'forces_{key}'])
            for key in keys
            if source.get(f'forces_{key}' is not None)
        ]

    def get_eigenvalues(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        eigenvalues = source.get('band_energies')
        if eigenvalues is None:
            return []
        n_spin = self.get_n_spin_channels()
        n_eigs = len(eigenvalues[0])
        n_bands = np.size(eigenvalues) // int(n_spin * n_eigs)
        eigenvalues = np.reshape(eigenvalues, (n_spin, n_bands, n_eigs)) * ureg.eV
        results = [dict(eigenvalues=eig) for n, eig in enumerate(eigenvalues)]
        occupations = source.get('occupation_numbers')
        if occupations is not None:
            occupations = np.reshape(occupations, (n_spin, n_bands, n_eigs))
            for n, occ in enumerate(occupations):
                results[n]['occupations'] = occ
        return results

    def get_configurations(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        methods = {
            'self_consistent': 'single_point',
            'bandstructure': 'single_point',
            'bfgs_geometry_optimization': 'geometry_optimization',
            'molecular_dynamics': 'molecular_dynamics',
            'langevin_dynamics': 'langevin_dynamics',
            'damped_dynamics': 'geometry_optimization',
            'vcs_wentzcovitch_damped_minimization': 'geometry_optimization',
        }

        configurations = []
        for key in methods:
            config = source.get(key)
            if config is None:
                continue
            configurations.append(config.get('self_consistent', config))
        return configurations


class PWSCFMainfileXMLParser(MainfileXMLParser):
    def get_configurations(self, source: dict[str, Any]):
        keys = ['input', 'output']
        return [source[key] for key in keys if source.get(key) is not None]


class PWSCFArchiveWriter(QuantumEspressoArchiveWriter):
    schema = pwscf
    _text_parser = PWSCFMainfileTextParser(text_parser=PWSCFFileParser())
    _xml_parser = PWSCFMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
