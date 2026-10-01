"""Cross-parser regression: FHI-aims `numerical_settings` must survive other
plugins' schema packages being loaded in the same process.

Per-parser tests import one schema module at a time and therefore never saw
the `mapper_m_def` poisoning: with the orca schema package imported, an
FHI-aims parse lost its typed top-level `KSpace`/`SelfConsistency` sections
(degraded to bare `NumericalSettings`), order-dependently. These tests run
the multi-plugin configurations CI previously never exercised.

Note on the legacy failure mode these tests guard against: the old reader
consumed the slot with `pop()`, so only the *first* parse in a process
exhibited the degradation — in a buggy world, whichever of these tests runs
first fails and shields the later one. One failure is all CI needs.
"""

import importlib

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

LOGGER = get_logger(__name__)
MAINFILE = 'tests/data/fhiaims/Si_geomopt/out.out'


def assert_numerical_settings_intact(archive):
    settings = archive.data.model_method[0].numerical_settings
    by_type = sorted(type(s).__name__ for s in settings)
    assert by_type == [
        'KSpace',
        'SelfConsistency',
        'SelfConsistency',
        'SelfConsistency',
    ], by_type

    kspace = next(s for s in settings if type(s).__name__ == 'KSpace')
    assert list(kspace.k_mesh[0].grid) == [8, 8, 8]
    for criterion in settings:
        if type(criterion).__name__ == 'SelfConsistency':
            assert criterion.threshold_change is not None


def test_fhiaims_numerical_settings_survive_orca_import():
    """The bisect repro as a test: orca's schema package annotates the shared
    `ModelMethod.numerical_settings` with `m_def=LocalCorrelationSettings`;
    importing it must not affect an FHI-aims parse. (Import order within one
    process is fixed by module caching, so the poisoning order — writer
    imported before the victim builds its mappers — is the one exercised.)"""
    importlib.import_module('nomad_simulation_parsers.schema_packages.orca')
    importlib.import_module('nomad_simulation_parsers.schema_packages.fhiaims')

    archive = EntryArchive()
    FHIAimsParser().parse(MAINFILE, archive, LOGGER)
    assert_numerical_settings_intact(archive)


def test_fhiaims_through_full_plugin_loading():
    """The production path: `nomad.client.parse` loads every plugin entry
    point, exactly like `nomad parse` or an Oasis."""
    from nomad.client import parse

    archive = parse(MAINFILE)[0]
    assert_numerical_settings_intact(archive)
