import functools
import inspect
import os
import re
from collections.abc import Callable
from glob import glob
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )
from nomad.utils import get_logger

DEFAULT_LOGGER = get_logger(__name__)


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
