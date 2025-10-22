RE_FLOAT = r'[-+]?\d+\.*\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'


def to_float(string: str | None) -> float | None:
    if string is None:
        return None
    try:
        value = float(string)
    except ValueError:
        value = None
    return value
