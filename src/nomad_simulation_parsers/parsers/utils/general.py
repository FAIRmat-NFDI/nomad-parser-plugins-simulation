import ast
import functools
import inspect
import os
import re
import textwrap
from collections.abc import Callable, Iterable
from glob import glob
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )
from nomad.utils import get_logger

DEFAULT_LOGGER = get_logger(__name__)

# Electronic structure constants
OCCUPATION_THRESHOLD = 0.5  # Threshold for occupied vs unoccupied states


def create_mapping_table(  # noqa: PLR0915
    file_parser: Any,
    archive_parser: Any,
    function_objects: Iterable[Any] | Any | None = None,
) -> list[dict[str, Any]]:
    """Compare file-parser quantities with the archive mapper's source paths.

    ``file_parser`` can be a ``TextParser`` or a parser wrapper exposing its
    ``text_parser``.  Nested ``Quantity`` objects are expanded using their
    ``sub_parser``.  ``archive_parser`` is the ``MetainfoParser`` used in the
    conversion.  Its mapper is inspected directly, including paths passed to
    transformation functions, so custom mappings such as
    ``get_eigenvalues(.eigenvalues_occupancies)`` are included. If
    ``function_objects`` is provided, functions named by the mapper are also
    inspected for direct accesses such as ``source.get('energy_total')``,
    ``source['energy_total']``, and ``source.energy_total``.

    Returns rows with ``quantity``, ``mapped``, and ``mapping`` keys.  The latter
    contains the mapper source paths that consume the file quantity.
    """

    parser = getattr(file_parser, 'text_parser', None) or file_parser
    quantities = getattr(parser, 'quantities', [])

    def iter_quantities(items: Iterable[Any], prefix: str = ''):
        for quantity in items:
            name = getattr(quantity, 'name', None)
            if not name:
                continue
            path = f'{prefix}.{name}' if prefix else name
            yield path, quantity
            sub_parser = getattr(quantity, 'sub_parser', None)
            if sub_parser is not None:
                yield from iter_quantities(getattr(sub_parser, 'quantities', []), path)

    def path_string(path: Any) -> str | None:
        if path is None:
            return None
        value = getattr(path, 'absolute_path', None) or getattr(path, 'path', None)
        return str(value) if value is not None else str(path)

    def mapper_paths(mapper: Any) -> list[str]:
        paths: list[str] = []
        paths.extend(getattr(mapper, 'all_paths', []))
        source = getattr(mapper, 'source', None)
        source_path = path_string(getattr(source, 'path', None))
        if source_path:
            paths.append(source_path)
        source_transformer = getattr(source, 'transformer', None)
        if source_transformer is not None:
            paths.extend(
                path
                for path in (
                    path_string(argument)
                    for argument in getattr(source_transformer, 'function_args', [])
                )
                if path
            )
        paths.extend(
            path
            for path in (
                path_string(argument)
                for argument in getattr(mapper, 'function_args', [])
            )
            if path
        )
        for child in getattr(mapper, 'mappers', []):
            paths.extend(mapper_paths(child))
        return paths

    def mapper_function_names(mapper: Any) -> set[str]:
        names = set()
        source_transformer = getattr(
            getattr(mapper, 'source', None), 'transformer', None
        )
        if source_transformer is not None and source_transformer.function_name:
            names.add(source_transformer.function_name)
        function_name = getattr(mapper, 'function_name', None)
        if function_name:
            names.add(function_name)
        for child in getattr(mapper, 'mappers', []):
            names.update(mapper_function_names(child))
        return names

    def function_accesses(objects: Iterable[Any] | Any | None) -> dict[str, set[str]]:
        if objects is None:
            return {}
        if isinstance(objects, type | str) or not isinstance(objects, Iterable):
            objects = [objects]
        objects = list(objects)
        accesses: dict[str, set[str]] = {}
        inspected: set[tuple[str, int]] = set()

        def find_function(name: str) -> Any:
            for obj in objects:
                candidate = getattr(obj, name, None)
                if callable(candidate):
                    return candidate
            return None

        def inspect_function(name: str, candidate: Any) -> set[str]:
            identity = (name, id(candidate))
            if identity in inspected:
                return set()
            inspected.add(identity)

            try:
                source = inspect.getsource(candidate)
            except (OSError, TypeError):
                return set()

            source = textwrap.dedent(source)
            try:
                function = ast.parse(source).body[0]
                arguments = function.args  # type: ignore[attr-defined]
                access_names = [
                    argument.arg
                    for argument in (
                        *getattr(arguments, 'posonlyargs', []),
                        *arguments.args,
                        *arguments.kwonlyargs,
                    )
                    if argument.arg not in {'self', 'cls'}
                ]
                access_names.extend(
                    node.id
                    for node in ast.walk(function)
                    if isinstance(node, ast.Name)
                    and isinstance(node.ctx, ast.Store)
                    and node.id not in {'self', 'cls'}
                )
            except (AttributeError, SyntaxError, TypeError):
                access_names = ['source']

            keys: set[str] = set()
            for access_name in set(access_names):
                parameter = re.escape(access_name)
                matches = re.findall(
                    rf"(?:\b{parameter}\s*\.\s*get\s*\(\s*['\"]([^'\"]+)"
                    rf"|\b{parameter}\s*\[\s*['\"]([^'\"]+)"
                    rf'|\b{parameter}\s*\.\s*(?!get\b)([A-Za-z_]\w*))',
                    source,
                )
                keys.update(key for match in matches for key in match if key)

            helper_names = {
                call.func.attr
                for call in ast.walk(function)
                if isinstance(call, ast.Call)
                and isinstance(call.func, ast.Attribute)
                and isinstance(call.func.value, ast.Name)
                and call.func.value.id in {'self', 'cls'}
            }
            helper_names.update(
                call.func.id
                for call in ast.walk(function)
                if isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
            )
            for helper_name in helper_names:
                helper = find_function(helper_name)
                if helper is not None:
                    keys.update(inspect_function(helper_name, helper))
            return keys

        for name in mapper_function_names(archive_parser.mapper):
            candidate = find_function(name)
            if candidate is None:
                continue
            accesses[name] = inspect_function(name, candidate)
        return accesses

    def normalize(path: str) -> str:
        path = path.strip().lstrip('.')
        path = path.replace('"', '').replace("'", '')
        path = re.sub(r'\[[^]]*\]', '', path)
        return path.strip('.')

    source_paths = [
        normalize(part)
        for path in mapper_paths(archive_parser.mapper)
        for part in path.split('||')
    ]
    source_paths = sorted(
        dict.fromkeys(path for path in source_paths if path and path != '@')
    )
    accessed_by = function_accesses(function_objects)
    function_names = sorted(mapper_function_names(archive_parser.mapper))
    table = []
    for quantity_path, _ in iter_quantities(quantities):
        normalized_quantity = normalize(quantity_path)
        consumed_by = [
            source
            for source in source_paths
            if source == normalized_quantity
            or source.startswith(f'{normalized_quantity}.')
            or normalized_quantity.startswith(f'{source}.')
            or normalized_quantity.endswith(f'.{source}')
        ]
        leaf = normalized_quantity.rsplit('.', 1)[-1]
        consumed_by.extend(
            f'{name}({leaf})'
            for name in function_names
            if leaf in accessed_by.get(name, set())
        )
        consumed_by = list(dict.fromkeys(consumed_by))
        table.append(
            {
                'quantity': quantity_path,
                'mapped': bool(consumed_by),
                'mapping': consumed_by,
            }
        )
    return table


def as_list(value: Any) -> list[Any]:
    """Wrap scalar payloads in a list; XML-derived sources are scalar for
    single occurrences and lists for repeated ones."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def link_outputs_to_model_systems(simulation) -> None:
    """Ensure output sections have stable step indices and model_system references.

    Strategy:
    - assign `output.step` when missing
    - assign `output.model_system_ref` when missing
      * index-matched if number of outputs equals number of model_systems
      * otherwise fallback to the last model_system when available
    """
    model_systems = simulation.model_system or []
    outputs = simulation.outputs or []

    for index, output in enumerate(outputs):
        if 'step' in output.m_def.all_quantities and output.step is None:
            output.step = index

        if output.model_system_ref is None:
            if len(model_systems) == len(outputs) and len(model_systems) > 0:
                output.model_system_ref = model_systems[index]
            elif len(model_systems) > 0:
                output.model_system_ref = model_systems[-1]


def calculate_band_gap_from_occupations(
    eigenvalues: np.ndarray | Any,
    occupations: np.ndarray | Any,
    occupation_threshold: float = OCCUPATION_THRESHOLD,
    spin_channel: int | None = None,
    energy_units: Any = None,
) -> dict[str, Any] | None:
    """Calculate band gap from eigenvalues and occupation numbers.

    This utility consolidates the common band gap calculation pattern used across
    multiple parsers (ABINIT, VASP, GPAW, Octopus, AMS, Exciting). The algorithm:
    1. Separates occupied from unoccupied states based on occupation threshold
    2. Finds highest occupied energy (valence band maximum)
    3. Finds lowest unoccupied energy (conduction band minimum)
    4. Computes gap = CBM - VBM (forced to be >= 0 for metals)

    Args:
        eigenvalues: Energy eigenvalues array (can have units via pint)
        occupations: Occupation numbers (typically 0-2 for spin-polarized,
            0-1 otherwise)
        occupation_threshold: Threshold to differentiate occupied/unoccupied
            (default: 0.5)
        spin_channel: Optional spin channel index (0=up, 1=down)
        energy_units: Optional pint units to apply to gap value (if
            eigenvalues are unitless)

    Returns:
        Dictionary with 'value' (and optional 'spin_channel') or None if
            invalid/insufficient data

    Example:
        >>> eigenvalues = np.array([[-5.0, -4.0, -3.0, 3.0, 4.0]])  # eV
        >>> occupations = np.array([[2.0, 2.0, 2.0, 0.0, 0.0]])
        >>> result = calculate_band_gap_from_occupations(eigenvalues, occupations)
        >>> result['value']  # Should be 6.0 (gap between -3.0 and 3.0)
        6.0
    """
    if eigenvalues is None or occupations is None:
        return None

    # Handle pint quantities by extracting values and units
    eigenvalues_values = eigenvalues
    eigenvalues_units = None
    if hasattr(eigenvalues, 'magnitude'):
        eigenvalues_values = eigenvalues.magnitude
        eigenvalues_units = eigenvalues.units

    # Convert to numpy arrays
    eigenvalues_arr = np.asarray(eigenvalues_values, dtype=float)
    occupations_arr = np.asarray(occupations, dtype=float)

    # Validate array properties
    if eigenvalues_arr.size == 0 or occupations_arr.size == 0:
        return None
    if eigenvalues_arr.shape != occupations_arr.shape:
        return None

    # Separate occupied and unoccupied states
    occupied = eigenvalues_arr[occupations_arr >= occupation_threshold]
    unoccupied = eigenvalues_arr[occupations_arr < occupation_threshold]

    # Need both occupied and unoccupied states to have a gap
    if occupied.size == 0 or unoccupied.size == 0:
        return None

    # Calculate gap (VBM to CBM)
    valence_max = np.max(occupied)
    conduction_min = np.min(unoccupied)
    gap = float(conduction_min - valence_max)

    # Force non-negative (metals have gap=0)
    gap = max(0.0, gap)

    # Apply units (prefer eigenvalues_units, fallback to energy_units parameter)
    if eigenvalues_units is not None:
        gap_value = gap * eigenvalues_units
    elif energy_units is not None:
        gap_value = gap * energy_units
    else:
        gap_value = gap

    # Build result dictionary
    result = {'value': gap_value}
    if spin_channel is not None:
        result['spin_channel'] = spin_channel

    return result


def search_files(pattern: str, basedir: str, **kwargs) -> list[str]:
    """Search files following the `pattern` starting from `basedir`. The search is
    performed recursively in all sub-folders (deep=True) or parent folders (deep=False).
    A futher regex search with `re_pattern` is done to filter the matching files.

    Args:
        pattern (str): pattern to match the files in the folder
        basedir (str): directory to start the search
        **deep (bool, optional): folders search direction (True=down, False=up)
        **re_pattern (str, optional): additional regex pattern to filter matching files
        **include_all (bool, optional): if True will include all matched files in sub
            directories to to max_dirs

    Returns:
        list: list of matching files
    """
    deep = kwargs.get('deep', True)
    max_dirs = kwargs.get('max_dirs', 10)
    re_pattern = kwargs.get('re_pattern', '')
    include_all = kwargs.get('include_all', False)

    filenames = []
    for _ in range(max_dirs):
        filenames.extend(glob(f'{basedir}/{pattern}'))
        pattern = os.path.join('**' if deep else '..', pattern)
        if filenames and not include_all:
            break

    if len(filenames) > 1 and re_pattern:
        # filter files that match
        matches = [f for f in filenames if re.search(re_pattern, f)]
        filenames = matches if matches else filenames

    filenames = [f for f in filenames if os.access(f, os.F_OK)]
    return filenames


def log(
    function: 'Callable' = None,
    logger: 'BoundLogger' = DEFAULT_LOGGER,
    exc_msg: str = None,
    exc_raise: bool = False,
    default: Any = None,
):
    """
    Function decorator to log exceptions.

    Args:
        function (Callable): function to evaluate
        logger (Logger, optional): logger to attach exceptions
        exc_msg (str, optional): prefix to exception
        exc_raise (bool, optional): if True will raise error
        default (Any, optional): return value of function if error
    """

    def _log(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            _logger = kwargs.get('logger', logger)
            _exc_msg = kwargs.get(
                'exc_msg', exc_msg or f'Exception raised in {func.__name__}:'
            )
            _exc_raise = kwargs.get('exc_raise', exc_raise)
            func.__annotations__['logger'] = _logger
            try:
                return func(
                    *args,
                    **{
                        key: val
                        for key, val in kwargs.items()
                        if key in inspect.signature(func).parameters
                    },
                )
            except Exception as e:
                _logger.warning(f'{_exc_msg} {e}')
                if _exc_raise:
                    raise e
                return kwargs.get('default', default)

        return wrapper

    return _log(function) if function else _log
