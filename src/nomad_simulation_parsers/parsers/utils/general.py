import functools
import inspect
import os
import re
from collections.abc import Callable
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
