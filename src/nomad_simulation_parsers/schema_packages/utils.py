from typing import Any

from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import Quantity, Section, SubSection
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY


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
