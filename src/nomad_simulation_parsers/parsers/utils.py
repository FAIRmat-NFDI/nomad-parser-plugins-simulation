import os
import re
from glob import glob


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
