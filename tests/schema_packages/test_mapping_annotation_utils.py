"""Unit tests for `add_mapping_annotation`'s `m_def=` handling.

The `m_def=` path must re-target the keyed annotation onto the subclass and
leave no trace on the shared definition: the removed
`Section.more['mapper_m_def']` write leaked across parser plugins (the
mapping parser consumed it for whichever parser built its mappers first,
stripping e.g. FHI-aims' `KSpace`/`SelfConsistency` whenever the orca schema
package was imported).
"""

import importlib

import pytest
from nomad_file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import general, model_method, numerical_settings

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

TEST_KEY = '_test_mdef'


@pytest.fixture(autouse=True)
def _cleanup_test_key():
    """Annotations persist on the shared nomad_simulations definitions for
    the lifetime of the process; strip the throwaway key after each test."""
    yield
    for definition in (
        model_method.ModelMethod.numerical_settings,
        numerical_settings.Smearing.m_def,
    ):
        definition.m_annotations.get(MAPPING_ANNOTATION_KEY, {}).pop(TEST_KEY, None)
    model_method.ModelMethod.numerical_settings.more.pop('mapper_m_def', None)


def test_m_def_retargets_onto_subclass():
    sub_section = model_method.ModelMethod.numerical_settings
    add_mapping_annotation(
        sub_section, TEST_KEY, '.@', m_def=numerical_settings.Smearing.m_def
    )

    smearing_annotations = numerical_settings.Smearing.m_def.m_annotations.get(
        MAPPING_ANNOTATION_KEY, {}
    )
    assert TEST_KEY in smearing_annotations
    assert smearing_annotations[TEST_KEY].mapper == '.@'
    # the subsection itself stays bare: resolution goes through the
    # polymorphic inheriting-section scan, scoped by the annotation key
    assert TEST_KEY not in sub_section.m_annotations.get(MAPPING_ANNOTATION_KEY, {})
    # deprecation guard: no unkeyed side channel on the shared definition
    assert 'mapper_m_def' not in sub_section.more


def test_m_def_fall_through_for_non_inheriting_section():
    sub_section = model_method.ModelMethod.numerical_settings
    # Program is not a subclass of NumericalSettings: the keyed annotation
    # falls through onto the subsection itself (documented behavior)
    add_mapping_annotation(sub_section, TEST_KEY, '.@', m_def=general.Program.m_def)

    assert TEST_KEY in sub_section.m_annotations.get(MAPPING_ANNOTATION_KEY, {})
    assert 'mapper_m_def' not in sub_section.more


def test_no_slot_after_importing_all_schema_packages():
    """Guards every current and future in-repo `m_def=` caller at once."""
    import nomad_simulation_parsers.schema_packages as schema_packages

    entry_points = [
        value
        for value in list(vars(schema_packages).values())
        if isinstance(value, schema_packages.EntryPoint)
    ]
    for entry_point in entry_points:
        importlib.import_module(entry_point.module)

    for definition in (
        model_method.ModelMethod.numerical_settings,
        general.Simulation.model_method,
        general.Simulation.outputs,
        general.Simulation.model_system,
    ):
        assert 'mapper_m_def' not in definition.more, definition
