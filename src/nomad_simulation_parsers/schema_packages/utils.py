from typing import Any, Iterable

from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import Quantity, Section, SubSection
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY


def remove_mapping_annotations(
    property: Section | SubSection,
    max_depth: int = 5,
    annotation_keys: Iterable[str] | None = None,
) -> None:
    """
    Remove mapping annotations from the input section definition, all its quantities
    and sub-sections recursively.

    Args:
        property: The section definition or subsection to remove annotations from.
        max_depth: The maximum depth of the recursion for sub-sections using the same
            section as parent.
        annotation_keys: Optional set of annotation keys that should be removed.
            When omitted, all mapping annotations are removed (legacy behaviour).
    """

    annotation_keys = set(annotation_keys) if annotation_keys is not None else None
    visited_sections: set[int] = set()

    def _clear_annotations(target: Section | SubSection | Quantity) -> None:
        mapping = target.m_annotations.get(MAPPING_ANNOTATION_KEY)
        if not mapping:
            return
        if annotation_keys is None:
            target.m_annotations.pop(MAPPING_ANNOTATION_KEY, None)
            return
        for key in annotation_keys:
            mapping.pop(key, None)
        if not mapping:
            target.m_annotations.pop(MAPPING_ANNOTATION_KEY, None)

    def _remove(property: Section | SubSection, depth: int = 0) -> None:
        if depth > max_depth:
            return

        _clear_annotations(property)

        property_section = (
            property.sub_section if isinstance(property, SubSection) else property
        )
        if property_section is None:
            return

        _clear_annotations(property_section)

        section_id = id(property_section)
        if section_id in visited_sections:
            return
        visited_sections.add(section_id)

        next_depth = depth + 1

        for quantity in property_section.all_quantities.values():
            _clear_annotations(quantity)

        for sub_section in property_section.all_sub_sections.values():
            _remove(sub_section, next_depth)

        for inheriting_section in property_section.all_inheriting_sections:
            _remove(inheriting_section, next_depth)

    _remove(property)


def add_mapping_annotation(
    property: Section | Quantity | SubSection,
    annotation_key: str,
    mapper: str | tuple[str, list[str] | tuple[str, list[str], dict[str, Any]]],
    update: bool = True,
    **kwargs,
) -> None:
    annotation = {annotation_key: Mapper(mapper=mapper, **kwargs)}
    if update:
        property.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(annotation)
    else:
        property.m_annotations[MAPPING_ANNOTATION_KEY] = annotation
