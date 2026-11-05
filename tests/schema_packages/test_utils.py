"""Unit tests for `add_mapping_annotation` and `add_mapping_annotations`.

The `m_def=` tests below guard a regression:
the `m_def=` path must re-target the keyed annotation onto the subclass and
leave no trace on the shared definition: the removed
`Section.more['mapper_m_def']` write leaked across parser plugins (the
mapping parser consumed it for whichever parser built its mappers first,
stripping e.g. FHI-aims' `KSpace`/`SelfConsistency` whenever the orca schema
package was imported).
"""

import importlib

import pytest
from nomad.datamodel import ArchiveSection
from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import Quantity
from nomad_file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import general, model_method, numerical_settings

from nomad_simulation_parsers import schema_packages
from nomad_simulation_parsers.schema_packages.utils import (
    add_mapping_annotation,
    add_mapping_annotations,
)

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


def make_target() -> type[ArchiveSection]:
    """A fresh throwaway section per call, so tests never annotate the shared
    `nomad_simulations` definitions and two calls never share state."""

    class Target(ArchiveSection):
        a = Quantity(type=str)
        b = Quantity(type=str)

    return Target


def mapping_annotations(target: type[ArchiveSection]) -> dict[str, dict]:
    definitions = {'m_def': target.m_def} | dict(target.m_def.all_quantities)
    return {
        name: dict(definition.m_annotations.get(MAPPING_ANNOTATION_KEY, {}))
        for name, definition in definitions.items()
    }


# The grouped call is defined as the single call applied per entry and per key,
# so the single call serves as the oracle and only the side effects (the stored
# annotations) are compared. Sketch for later, not adopted in this repo yet: a
# Hypothesis strategy drawing keys, mappers and kwargs could replace these
# hand-picked cases.
@pytest.mark.parametrize(
    'keys, mappers, kwargs',
    [
        ('k1', {'a': '.a'}, {}),
        (['k1', 'k2'], {'a': '.a', 'b': ('get_b', ['.b'])}, {}),
        (('k1', 'k2'), {'a': '.a'}, {}),
        (
            ['k1', 'k2'],
            {'a': '.a', 'b': '.b'},
            {'update_mode': 'append', 'cache': True},
        ),
        ('k1', {'m_def': '.@', 'a': '.a'}, {}),
        ('k1', {}, {}),
    ],
    ids=[
        'single-key',
        'key-list-with-transformer',
        'key-tuple',
        'shared-kwargs',
        'section-and-quantity',
        'empty-mappers',
    ],
)
def test_grouped_call_equals_single_calls(keys, mappers, kwargs):
    grouped, looped = make_target(), make_target()
    for target in (grouped, looped):
        add_mapping_annotation(target.a, 'k0', '.existing')

    add_mapping_annotations(
        keys, {getattr(grouped, name): m for name, m in mappers.items()}, **kwargs
    )
    for name, mapper in mappers.items():
        for key in [keys] if isinstance(keys, str) else keys:
            add_mapping_annotation(getattr(looped, name), key, mapper, **kwargs)

    assert mapping_annotations(grouped) == mapping_annotations(looped)


def test_key_list_gets_one_mapper_per_key():
    target = make_target()
    add_mapping_annotation(target.a, ['k1', 'k2'], '.a')

    stored = target.a.m_annotations[MAPPING_ANNOTATION_KEY]
    assert stored['k1'] is not stored['k2']


def test_key_list_with_update_false_replaces_and_keeps_all_keys():
    target = make_target()
    add_mapping_annotation(target.a, 'k0', '.existing')
    add_mapping_annotation(target.a, ['k1', 'k2'], '.a', update=False)

    stored = target.a.m_annotations[MAPPING_ANNOTATION_KEY]
    assert sorted(stored) == ['k1', 'k2']


def test_empty_key_list_raises():
    with pytest.raises(ValueError):
        add_mapping_annotation(make_target().a, [], '.a')


@pytest.mark.parametrize(
    'entry, kwargs, expected',
    [
        (
            Mapper(mapper='.a', update_mode='replace'),
            {'update_mode': 'append'},
            {'mapper': '.a', 'update_mode': 'replace', 'cache': None},
        ),
        (
            Mapper(mapper='.a'),
            {'update_mode': 'append', 'cache': True},
            {'mapper': '.a', 'update_mode': 'append', 'cache': True},
        ),
    ],
    ids=['entry-field-overrides-shared', 'unset-entry-fields-keep-shared'],
)
def test_mapper_entry_precedence(entry, kwargs, expected):
    target = make_target()
    add_mapping_annotations('k1', {target.a: entry}, **kwargs)

    stored = target.a.m_annotations[MAPPING_ANNOTATION_KEY]['k1']
    assert stored.model_dump(include=set(expected)) == expected
