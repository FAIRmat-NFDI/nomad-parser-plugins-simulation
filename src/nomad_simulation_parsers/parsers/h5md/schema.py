from __future__ import annotations

from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np

try:
    import h5py
except ImportError:
    h5py = None

try:
    import yaml
except ImportError:
    yaml = None


DEFAULT_SCHEMA = 'h5md-nomad.schema.yaml'


@dataclass
class ValidationIssue:
    path: str
    message: str


@dataclass
class ValidationResult:
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return not self.issues

    def add(self, path: str, message: str) -> None:
        self.issues.append(ValidationIssue(path=path, message=message))


def load_schema(path: str | Path | None = None) -> dict[str, Any]:
    """Load a H5MD schema DSL document from YAML."""
    if yaml is None:
        raise RuntimeError('PyYAML is required to load H5MD schema files.')

    if path is None:
        schema_text = (
            resources.files('nomad_simulation_parsers.parsers.h5md.schemas')
            .joinpath(DEFAULT_SCHEMA)
            .read_text()
        )
        return yaml.safe_load(schema_text)

    with Path(path).open() as schema_file:
        return yaml.safe_load(schema_file)


def validate_hdf5_file(
    path: str | Path, schema: dict[str, Any] | None = None
) -> ValidationResult:
    """Validate an HDF5 file against a H5MD schema DSL document."""
    if h5py is None:
        raise RuntimeError('h5py is required to validate HDF5 files.')

    schema = schema or load_schema()
    with h5py.File(path, 'r') as h5_file:
        return validate_hdf5(h5_file, schema)


def validate_hdf5(h5_object: Any, schema: dict[str, Any]) -> ValidationResult:
    """Validate an open h5py file or group against a H5MD schema DSL document."""
    if h5py is None:
        raise RuntimeError('h5py is required to validate HDF5 files.')

    result = ValidationResult()
    dimensions: dict[str, int] = {}
    _validate_node(h5_object, schema['root'], '/', result, dimensions)
    return result


def _validate_node(
    h5_object: Any,
    node_schema: dict[str, Any],
    path: str,
    result: ValidationResult,
    dimensions: dict[str, int],
) -> None:
    node_type = node_schema.get('type', 'group')
    if node_type == 'group':
        _validate_group(h5_object, node_schema, path, result, dimensions)
    elif node_type == 'dataset':
        _validate_dataset(h5_object, node_schema, path, result, dimensions)
    else:
        result.add(path, f'Unsupported node type {node_type!r}.')


def _validate_group(
    group: Any,
    group_schema: dict[str, Any],
    path: str,
    result: ValidationResult,
    dimensions: dict[str, int],
) -> None:
    if not isinstance(group, h5py.Group):
        result.add(path, 'Expected HDF5 group.')
        return

    _validate_attributes(
        group.attrs,
        group_schema.get('attributes', {}),
        path,
        result,
        dimensions,
    )

    for name, child_schema in group_schema.get('children', {}).items():
        child_path = _join_path(path, name)
        required = child_schema.get('required', True)
        child_type = child_schema.get('type', 'group')

        if child_type in {'soft_link', 'external_link'}:
            _validate_link(group, name, child_schema, child_path, result)
            continue

        if name not in group:
            if required:
                result.add(child_path, 'Missing required child.')
            continue

        _validate_node(group[name], child_schema, child_path, result, dimensions)


def _validate_dataset(
    dataset: Any,
    dataset_schema: dict[str, Any],
    path: str,
    result: ValidationResult,
    dimensions: dict[str, int],
) -> None:
    if not isinstance(dataset, h5py.Dataset):
        result.add(path, 'Expected HDF5 dataset.')
        return

    _validate_dtype(dataset.dtype, dataset_schema.get('dtype'), path, result)
    _validate_shape(
        dataset.shape, dataset_schema.get('shape'), path, result, dimensions
    )
    _validate_attributes(
        dataset.attrs,
        dataset_schema.get('attributes', {}),
        path,
        result,
        dimensions,
    )


def _validate_link(
    group: Any,
    name: str,
    link_schema: dict[str, Any],
    path: str,
    result: ValidationResult,
) -> None:
    required = link_schema.get('required', True)
    link = group.get(name, getlink=True)
    if link is None:
        if required:
            result.add(path, 'Missing required link.')
        return

    link_type = link_schema.get('type')
    expected_link_type = {
        'soft_link': h5py.SoftLink,
        'external_link': h5py.ExternalLink,
    }[link_type]
    if not isinstance(link, expected_link_type):
        result.add(path, f'Expected {link_type}.')


def _validate_attributes(
    attributes: Any,
    attribute_schemas: dict[str, Any],
    owner_path: str,
    result: ValidationResult,
    dimensions: dict[str, int],
) -> None:
    for name, attribute_schema in attribute_schemas.items():
        if owner_path == '/':
            attribute_path = f'{owner_path}@{name}'
        else:
            attribute_path = f'{owner_path}/@{name}'
        required = attribute_schema.get('required', True)
        if name not in attributes:
            if required:
                result.add(attribute_path, 'Missing required attribute.')
            continue

        value = attributes[name]
        _validate_dtype(
            np.asarray(value).dtype,
            attribute_schema.get('dtype'),
            attribute_path,
            result,
        )
        _validate_shape(
            np.shape(value),
            attribute_schema.get('shape'),
            attribute_path,
            result,
            dimensions,
        )


def _validate_dtype(
    actual_dtype: Any,
    expected_dtype: str | None,
    path: str,
    result: ValidationResult,
) -> None:
    if expected_dtype in {None, 'any'}:
        return

    kind = np.dtype(actual_dtype).kind
    expected_kinds = {
        'bool': {'b'},
        'boolean': {'b'},
        'float': {'f'},
        'integer': {'i', 'u'},
        'numeric': {'f', 'i', 'u'},
        'string': {'O', 'S', 'U'},
    }

    if expected_dtype not in expected_kinds:
        result.add(path, f'Unsupported dtype constraint {expected_dtype!r}.')
        return

    if kind not in expected_kinds[expected_dtype]:
        result.add(path, f'Expected dtype {expected_dtype}, got {actual_dtype}.')


def _validate_shape(
    actual_shape: tuple[int, ...],
    expected_shape: list[Any] | None,
    path: str,
    result: ValidationResult,
    dimensions: dict[str, int],
) -> None:
    if expected_shape is None:
        return

    expected_shape = list(expected_shape)
    if len(actual_shape) != len(expected_shape):
        result.add(
            path,
            f'Expected rank {len(expected_shape)}, got rank {len(actual_shape)}.',
        )
        return

    for index, (actual, expected) in enumerate(zip(actual_shape, expected_shape)):
        if expected in {'*', 'any'}:
            continue
        if isinstance(expected, int):
            if actual != expected:
                result.add(path, f'Expected shape[{index}]={expected}, got {actual}.')
            continue
        if isinstance(expected, str):
            previous = dimensions.setdefault(expected, actual)
            if previous != actual:
                result.add(
                    path,
                    f'Expected dimension {expected!r}={previous}, got {actual}.',
                )
            continue
        result.add(path, f'Unsupported shape constraint {expected!r}.')


def _join_path(parent: str, child: str) -> str:
    if parent == '/':
        return f'/{child}'
    return f'{parent}/{child}'
