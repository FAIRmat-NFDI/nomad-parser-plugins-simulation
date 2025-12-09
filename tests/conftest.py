from collections.abc import Sequence
from typing import Any

import magic
import numpy as np
import pint
import pytest
from nomad.config import config
from nomad.parsing.parsers import _compressions, encoding_magic


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


def compare_values(
    value: Any, reference: Any, all_keys: bool = False, check_type: bool = False
) -> bool:
    """
    Compare a value to reference recursively.
    Parameters:
        value[Any]: The value to compare.
        reference[Any]: The reference value to compare to.
        all_keys[bool]: For dictionaries, value and reference must have the same keys.
    """

    def compare_dicts(value: dict[str, Any], reference: dict[str, Any]) -> bool:
        equal = True
        if all_keys and sorted(list(value.keys())) != sorted(list(reference.keys())):
            return False

        for key, val in reference.items():
            equal = compare_values(value.get(key), val)
            if not equal:
                break
        return equal

    def compare_lists(value: list[Any], reference: list[Any]) -> bool:
        equal = True
        if len(reference) != len(value):
            return False

        types = (dict, np.ndarray, Sequence)
        if isinstance(reference[0], types) or isinstance(value[0], types):
            for n, val in enumerate(reference):
                equal = compare_values(value[n], val)
                if not equal:
                    break
        else:
            equal = value == reference
        return equal

    equal = True
    if isinstance(value, np.ndarray):
        equal = np.allclose(value, reference)

    elif isinstance(value, pint.Quantity):
        equal = compare_values(value.magnitude, reference.get('__value')) and str(
            value.units
        ) == reference.get('__unit')

    elif isinstance(value, float):
        equal = value == approx(reference)

    elif isinstance(reference, dict):
        equal = compare_dicts(value, reference)

    elif isinstance(reference, list) and reference:
        equal = compare_lists(value, reference)

    else:
        equal = value == reference
        if check_type and equal:
            equal = isinstance(value, type(reference))

    return equal


def get_child_archive_keys(filename, parser):
    with open(filename, 'rb') as f:
        compression, open_compressed = _compressions.get(f.read(3), (None, open))

    with open_compressed(filename, 'rb') as cf:  # type: ignore
        buffer = cf.read(config.process.parser_matching_size)

    mime_type = magic.from_buffer(buffer, mime=True)

    decoded_buffer = None
    encoding = None
    try:  # Try to open the file as a string for regex matching.
        decoded_buffer = buffer.decode('utf-8')
    except UnicodeDecodeError:
        # This file is either binary or has wrong encoding
        encoding = encoding_magic.from_buffer(buffer)

        if config.services.force_raw_file_decoding:
            encoding = 'iso-8859-1'

        if encoding in ['iso-8859-1']:
            try:
                decoded_buffer = buffer.decode(encoding)
            except Exception:
                pass

    return parser.is_mainfile(filename, mime_type, buffer, decoded_buffer, compression)
