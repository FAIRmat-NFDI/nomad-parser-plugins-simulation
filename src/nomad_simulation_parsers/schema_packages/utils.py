from typing import Any

from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import Quantity, Section, SubSection
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY


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
    annotation_key: str,
    mapper: str | tuple[str, list[str] | tuple[str, list[str], dict[str, Any]]],
    update: bool = True,
    m_def: Section | None = None,
    **kwargs,
) -> None:
    annotation = {annotation_key: Mapper(mapper=mapper, **kwargs)}
    if m_def is not None and isinstance(property, SubSection):
        property.more['mapper_m_def'] = m_def.qualified_name()
        for inheriting_section in property.sub_section.all_inheriting_sections or []:
            if m_def.qualified_name() == inheriting_section.qualified_name():
                add_mapping_annotation(inheriting_section, annotation_key, mapper, update, **kwargs)
                return

    if update:
        property.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(annotation)
    else:
        property.m_annotations[MAPPING_ANNOTATION_KEY] = annotation
