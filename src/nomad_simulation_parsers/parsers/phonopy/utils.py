import re
from pathlib import Path


def find_phonopy_yaml(mainfile: str, calculator: str) -> Path | None:
    """Find a Phonopy result beside or above the first displaced calculation."""
    directory = Path(mainfile).parent
    yaml_file = directory / 'phonopy.yaml'
    if not yaml_file.is_file() and re.fullmatch(
        r'(?:disp(?:lacement)?|supercell|phonopy-.*-displacement)[-_]0*1',
        directory.name,
        re.IGNORECASE,
    ):
        yaml_file = directory.parent / 'phonopy.yaml'
    if not yaml_file.is_file():
        return None
    try:
        with yaml_file.open(encoding='utf-8') as stream:
            header = stream.read(4096)
    except OSError:
        return None
    match = re.search(r'^  calculator:\s*[\'"]?([\w-]+)', header, re.MULTILINE)
    # Older VASP Phonopy files omit the calculator field.
    yaml_calculator = match.group(1).lower() if match else 'vasp'
    return yaml_file if yaml_calculator == calculator.lower() else None
