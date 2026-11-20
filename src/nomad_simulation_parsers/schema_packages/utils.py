from collections.abc import Sequence
from typing import Any

from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import Quantity, Section, SubSection
from nomad_file_parser.mapping_parser import MAPPING_ANNOTATION_KEY


def remove_mapping_annotations(property: Section, max_depth: int = 5) -> None:
    """
    Remove mapping annotations from the input section definition, all its quantities
    and sub-sections recursively.
    Args:
        property (Section): The section definition to remove the annotations from.
        max_depth (int, optional): The maximum depth of the recursion for sub-sections
            using the same section as parent.
    """

    def _remove(property: Section | SubSection, depth: int = 0):
        if depth > max_depth:
            return

        annotation_key = 'mapping'
        property.m_annotations.pop(annotation_key, None)

        depth += 1
        property_section = (
            property.sub_section if isinstance(property, SubSection) else property
        )
        for quantity in property_section.all_quantities.values():
            quantity.m_annotations.pop(annotation_key, None)

        for sub_section in property_section.all_sub_sections.values():
            if sub_section.m_annotations.get(annotation_key):
                _remove(sub_section, depth)
            elif sub_section.sub_section.m_annotations.get(annotation_key):
                _remove(sub_section.sub_section, depth)
            else:
                for (
                    inheriting_section
                ) in sub_section.sub_section.all_inheriting_sections:
                    if inheriting_section.m_annotations.get(annotation_key):
                        _remove(inheriting_section, depth)

    _remove(property)


def add_mapping_annotation(
    property: Section | Quantity | SubSection,
    annotation_key: str | Sequence[str],
    mapper: str | tuple[str, list[str] | tuple[str, list[str], dict[str, Any]]],
    update: bool = True,
    m_def: Section | None = None,
    **kwargs,
) -> None:
    """
    Bind a schema ``property`` to a source path/transformer under
    ``annotation_key``, consumed by the ``nomad_file_parser`` mapping parser to
    populate the archive from parsed data.

    ``annotation_key`` may also be a sequence of keys, which registers the same
    ``mapper`` (and the same ``m_def`` and ``kwargs``) under each of them, e.g. when
    two files of a code share a layout. Each key gets its own ``Mapper`` instance.
    With ``update=False``, the keys of the sequence replace all mapping annotations
    previously stored on ``property``, and they are kept together. An empty sequence
    raises a ``ValueError``.

    ``m_def`` only redirects where the annotation is stored. If ``property`` is a
    ``SubSection`` and ``m_def`` is one of the sections inheriting from its declared
    type, the annotation is stored on that section definition instead, as if the
    helper had been called with ``m_def`` as ``property``; the subsection itself is
    left without an annotation, and the mapping parser finds the subclass through
    its scan of inheriting sections under the same ``annotation_key``. If ``m_def``
    matches none of them, or ``property`` is not a ``SubSection``, ``m_def`` is
    ignored and the annotation is stored on ``property``.

    Caveat for polymorphic subsections: the mapping conversion builds subsections
    by instantiating their *declared* type
    (``section.sub_section.section_cls()`` in ``nomad_file_parser``
    ``file_parser.parse_section``). A concrete subclass instance assigned to a
    subsection before conversion is therefore rebuilt as the declared base type
    (e.g. ``EnergyConvergenceTarget``/``ForceConvergenceTarget`` collapse to
    ``WorkflowConvergenceTarget``). To keep a polymorphic subclass, populate it
    manually *after* the parser's ``convert()`` pass instead of via this
    annotation (see ``abinit``/``exciting`` ``parse_workflow``).
    """
    if not isinstance(annotation_key, str):
        if not annotation_key:
            raise ValueError('annotation_key must contain at least one key')
        for index, key in enumerate(annotation_key):
            add_mapping_annotation(
                property, key, mapper, update or index > 0, m_def, **kwargs
            )
        return

    annotation = {annotation_key: Mapper(mapper=mapper, **kwargs)}
    if m_def is not None and isinstance(property, SubSection):
        for inheriting_section in property.sub_section.all_inheriting_sections or []:
            if m_def.qualified_name() == inheriting_section.qualified_name():
                add_mapping_annotation(
                    inheriting_section, annotation_key, mapper, update, **kwargs
                )
                return

    if update:
        property.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(annotation)
    else:
        property.m_annotations[MAPPING_ANNOTATION_KEY] = annotation


def add_mapping_annotations(
    annotation_key: str | Sequence[str],
    mappers: dict[
        Section | Quantity | SubSection,
        str | tuple[str, list[str] | tuple[str, list[str], dict[str, Any]]] | Mapper,
    ],
    **kwargs,
) -> None:
    """
    Group the mapping annotations of one section under ``annotation_key``, in place
    of a wrapper class around single ``add_mapping_annotation`` calls.

    Each entry of ``mappers`` binds a definition (a quantity, a subsection, or a
    section via ``m_def``) to a mapper, and is registered with
    ``add_mapping_annotation``. The targets are definitions, not names, so a
    misspelled target fails at import. ``kwargs`` (e.g. ``update_mode``) apply to
    every entry; to set them for one entry only, give that entry a ``Mapper``, whose
    explicitly set fields take precedence over ``kwargs``.

    Example:
        add_mapping_annotations(
            XML_KEY,
            {
                general.Program.name: '.i[?"@name"==\\'program\\'] | [0].__value',
                general.Program.version: '.i[?"@name"==\\'version\\'] | [0].__value',
            },
        )
    """
    for property, mapper in mappers.items():
        if isinstance(mapper, Mapper):
            entry_kwargs = mapper.model_dump(
                include=mapper.model_fields_set - {'mapper'}
            )
            add_mapping_annotation(
                property, annotation_key, mapper.mapper, **(kwargs | entry_kwargs)
            )
        else:
            add_mapping_annotation(property, annotation_key, mapper, **kwargs)
