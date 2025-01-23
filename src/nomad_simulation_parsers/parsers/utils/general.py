import functools
import os
import re
from glob import glob
from typing import TYPE_CHECKING, Any, Callable, Union

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )
from nomad.metainfo import Section, SubSection
from nomad.utils import get_logger


def search_files(
    pattern: str,
    basedir: str,
    deep: bool = True,
    max_dirs: int = 10,
    re_pattern: str = '',
) -> list[str]:
    """Search files following the `pattern` starting from `basedir`. The search is
    performed recursively in all sub-folders (deep=True) or parent folders (deep=False).
    A futher regex search with `re_pattern` is done to filter the matching files.

    Args:
        pattern (str): pattern to match the files in the folder
        basedir (str): directory to start the search
        deep (bool, optional): folders search direction (True=down, False=up)
        re_pattern (str, optional): additional regex pattern to filter matching files

    Returns:
        list: list of matching files
    """

    for _ in range(max_dirs):
        filenames = glob(f'{basedir}/{pattern}')
        pattern = os.path.join('**' if deep else '..', pattern)
        if filenames:
            break

    if len(filenames) > 1:
        # filter files that match
        matches = [f for f in filenames if re.search(re_pattern, f)]
        filenames = matches if matches else filenames

    filenames = [f for f in filenames if os.access(f, os.F_OK)]
    return filenames


def remove_mapping_annotations(property: Section, max_depth: int = 5) -> None:
    """
    Remove mapping annotations from the input section definition, all its quantities
    and sub-sections recursively.

    Args:
        property (Section): The section definition to remove the annotations from.
        max_depth (int, optional): The maximum depth of the recursion for sub-sections
            using the same section as parent.
    """

    def _remove(property: Union[Section, SubSection], depth: int = 0):
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


def log(
    function: Callable = None,
    logger: 'BoundLogger' = get_logger(__name__),
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

            try:
                return func(*args, **kwargs)
            except Exception as e:
                _logger.warning(f'{_exc_msg} {e}')
                if _exc_raise:
                    raise e
                return kwargs.get('default', default)

        return wrapper

    return _log(function) if function else _log
