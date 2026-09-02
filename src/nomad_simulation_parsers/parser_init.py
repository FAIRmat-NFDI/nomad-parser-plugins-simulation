from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import yaml


def _to_snake_case(value: str) -> str:
    return value.strip().replace('-', '_').replace(' ', '_').lower()


def _to_pascal_case(value: str) -> str:
    return ''.join(
        part.capitalize() for part in _to_snake_case(value).split('_') if part
    )


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open('r', encoding='utf-8') as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError('Metadata YAML must contain a top-level mapping.')
    return data


def _render_file_parser_template(parser_class_name: str) -> str:
    return f"""from nomad_file_parser import Quantity, TextParser


class {parser_class_name}OutParser(TextParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self):
        self._quantities = [
            # TODO: Define Quantity extractors for the mainfile.
            # Quantity('key', r'regex', repeats=False),
        ]
"""


def _render_parser_template(
    parser_class_name: str, code_name: str, use_mapping_parser: bool
) -> str:
    if use_mapping_parser:
        return f"""from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from .file_parser import {parser_class_name}OutParser

LOGGER = get_logger(__name__)


class {parser_class_name}MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def __init__(self):
        super().__init__(text_parser={parser_class_name}OutParser())


class {parser_class_name}MetainfoParser(MetainfoParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class {parser_class_name}ArchiveWriter(ArchiveWriter):
    code_name = '{code_name}'
    mainfile_parser = {parser_class_name}MainfileParser()
    metainfo_parser = {parser_class_name}MetainfoParser()

    def write_to_archive(self):
        self.archive.data = Simulation(program=Program(name=self.code_name))
        self.metainfo_parser.data_object = self.archive.data
        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.convert(self.metainfo_parser)
        # TODO: Implement the archive writing logic.


class {parser_class_name}(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger) -> None:
        archive_writer = {parser_class_name}ArchiveWriter()
        archive_writer.write(mainfile, archive, logger)
"""

    return f"""from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad_file_parser import ArchiveWriter
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from .file_parser import {parser_class_name}OutParser


class {parser_class_name}ArchiveWriter(ArchiveWriter):
    code_name = '{code_name}'
    mainfile_parser = {parser_class_name}OutParser()

    def write_to_archive(self):
        self.archive.data = Simulation(program=Program(name=self.code_name))
        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.parse()
        # TODO: Map parsed mainfile data from self.mainfile_parser to self.archive.data.


class {parser_class_name}(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger) -> None:
        archive_writer = {parser_class_name}ArchiveWriter()
        archive_writer.write(mainfile, archive, logger)
"""


def _render_schema_template(out_key: str, use_mapping_parser: bool) -> str:
    if use_mapping_parser:
        return f"""from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general
from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = '{out_key}'

# TODO: Add class-level mapping annotations here.

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')

m_package.__init_metainfo__()
"""

    return """from nomad.metainfo import SchemaPackage

m_package = SchemaPackage()

m_package.__init_metainfo__()
"""


def _render_parser_entry_block(  # noqa: PLR0913
    *,
    entry_var_name: str,
    code_key: str,
    parser_description: str,
    parser_class_path: str,
    aliases: list[str],
    code_name: str,
    code_homepage: str | None,
    parser_options: dict[str, Any],
) -> str:
    lines = [
        f'{entry_var_name} = EntryPoint(',
        f"    name='parsers/{code_key}',",
        f'    aliases={aliases!r},',
        f'    description={parser_description!r},',
        "    python_package='nomad_simulation_parsers',",
        f"    parser_class_name='{parser_class_path}',",
        f'    code_name={code_name!r},',
    ]
    if code_homepage:
        lines.append(f'    code_homepage={code_homepage!r},')

    for key in (
        'mainfile_contents_re',
        'mainfile_name_re',
        'mainfile_mime_re',
        'mainfile_binary_header_re',
        'mainfile_contents_dict',
        'mainfile_alternative',
        'supported_compressions',
        'code_category',
    ):
        if key in parser_options:
            lines.append(f'    {key}={parser_options[key]!r},')

    lines.append(')')
    return '\n'.join(lines)


def _render_schema_entry_block(
    *,
    entry_var_name: str,
    schema_entry_name: str,
    schema_description: str,
    schema_module: str,
) -> str:
    return '\n'.join(
        [
            f'{entry_var_name} = EntryPoint(',
            f'    name={schema_entry_name!r},',
            f'    description={schema_description!r},',
            f'    module={schema_module!r},',
            ')',
        ]
    )


def _append_if_missing(file_path: Path, anchor_text: str, content: str) -> bool:
    text = file_path.read_text(encoding='utf-8')
    if anchor_text in text:
        return False
    if not text.endswith('\n'):
        text += '\n'
    text = f'{text}\n{content}\n'
    file_path.write_text(text, encoding='utf-8')
    return True


def _insert_pyproject_entry_lines(
    pyproject_path: Path,
    parser_key: str,
    parser_value: str,
    schema_key: str,
    schema_value: str,
) -> list[str]:
    text = pyproject_path.read_text(encoding='utf-8')
    section_header = "[project.entry-points.'nomad.plugin']"
    section_start = text.find(section_header)
    if section_start < 0:
        msg = "Could not find [project.entry-points.'nomad.plugin'] section."
        raise ValueError(msg)

    section_body_start = section_start + len(section_header)
    section_tail = text[section_body_start:]

    next_section_index = section_tail.find('\n[')
    if next_section_index == -1:
        section_end = len(text)
    else:
        section_end = section_body_start + next_section_index

    section_content = text[section_body_start:section_end]
    updated = False
    added: list[str] = []

    parser_line = f'{parser_key} = "{parser_value}"'
    schema_line = f'{schema_key} = "{schema_value}"'

    if parser_line not in section_content:
        section_content = f'{section_content.rstrip()}\n{parser_line}\n'
        updated = True
        added.append(parser_line)

    if schema_line not in section_content:
        section_content = f'{section_content.rstrip()}\n{schema_line}\n'
        updated = True
        added.append(schema_line)

    if updated:
        new_text = f'{text[:section_body_start]}{section_content}{text[section_end:]}'
        pyproject_path.write_text(new_text, encoding='utf-8')

    return added


def initialize_from_metadata(  # noqa: PLR0915
    metadata_path: Path,
    root: Path,
    overwrite: bool = False,
    dry_run: bool = False,
    use_mapping_parser: bool = True,
) -> list[str]:
    metadata = _read_yaml(metadata_path)

    code_key_raw = metadata.get('code_key') or metadata.get('name')
    if not code_key_raw:
        raise ValueError('Metadata requires either "code_key" or "name".')
    code_key = _to_snake_case(str(code_key_raw))

    code_name = str(metadata.get('name') or code_key_raw)
    code_homepage = metadata.get('code_homepage')
    parser_description = str(
        metadata.get('description') or f'NOMAD parser for {code_name}.'
    )

    parser_meta = (
        metadata.get('parser') if isinstance(metadata.get('parser'), dict) else {}
    )
    schema_meta = (
        metadata.get('schema_package')
        if isinstance(metadata.get('schema_package'), dict)
        else {}
    )

    parser_class_name = str(
        parser_meta.get('class_name') or f'{_to_pascal_case(code_key)}Parser'
    )
    parser_entry_var_name = str(parser_meta.get('entry_point') or f'{code_key}_parser')
    parser_aliases = parser_meta.get('aliases')
    if parser_aliases is None:
        parser_aliases = [f'parsers/{code_key}']
    if not isinstance(parser_aliases, list) or not all(
        isinstance(alias, str) for alias in parser_aliases
    ):
        raise ValueError('parser.aliases must be a list of strings.')

    parser_class_path = (
        f'nomad_simulation_parsers.parsers.{code_key}.parser.{parser_class_name}'
    )
    schema_entry_var_name = str(
        schema_meta.get('entry_point') or f'{code_key}_schema_package'
    )
    schema_entry_name = str(
        schema_meta.get('name') or f'{_to_pascal_case(code_key)}SchemaPackage'
    )
    schema_description = str(
        schema_meta.get('description') or f'Schema package for {code_name}.'
    )

    schema_module_suffix = str(schema_meta.get('module_suffix') or code_key)
    schema_module = f'nomad_simulation_parsers.schema_packages.{schema_module_suffix}'

    parser_options = {
        key: parser_meta[key]
        for key in (
            'mainfile_contents_re',
            'mainfile_name_re',
            'mainfile_mime_re',
            'mainfile_binary_header_re',
            'mainfile_contents_dict',
            'mainfile_alternative',
            'supported_compressions',
            'code_category',
        )
        if key in parser_meta
    }

    parser_dir = root / 'src' / 'nomad_simulation_parsers' / 'parsers' / code_key
    parser_file = parser_dir / 'parser.py'
    file_parser_file = parser_dir / 'file_parser.py'
    parser_init = parser_dir / '__init__.py'
    schema_file = (
        root
        / 'src'
        / 'nomad_simulation_parsers'
        / 'schema_packages'
        / f'{schema_module_suffix}.py'
    )
    parsers_registry = (
        root / 'src' / 'nomad_simulation_parsers' / 'parsers' / '__init__.py'
    )
    schema_registry = (
        root / 'src' / 'nomad_simulation_parsers' / 'schema_packages' / '__init__.py'
    )
    pyproject_path = root / 'pyproject.toml'

    changed: list[str] = []

    if parser_file.exists() and not overwrite:
        raise FileExistsError(
            f'{parser_file} already exists. Re-run with --overwrite to replace it.'
        )

    if file_parser_file.exists() and not overwrite:
        raise FileExistsError(
            f'{file_parser_file} already exists. Re-run with --overwrite to replace it.'
        )

    if schema_file.exists() and not overwrite:
        raise FileExistsError(
            f'{schema_file} already exists. Re-run with --overwrite to replace it.'
        )

    parser_content = _render_parser_template(
        parser_class_name, code_name, use_mapping_parser
    )
    file_parser_content = _render_file_parser_template(parser_class_name)
    schema_content = _render_schema_template(code_key, use_mapping_parser)
    parser_entry_block = _render_parser_entry_block(
        entry_var_name=parser_entry_var_name,
        code_key=code_key,
        parser_description=parser_description,
        parser_class_path=parser_class_path,
        aliases=parser_aliases,
        code_name=code_name,
        code_homepage=code_homepage,
        parser_options=parser_options,
    )
    schema_entry_block = _render_schema_entry_block(
        entry_var_name=schema_entry_var_name,
        schema_entry_name=schema_entry_name,
        schema_description=schema_description,
        schema_module=schema_module,
    )

    parser_pyproject_value = f'nomad_simulation_parsers.parsers:{parser_entry_var_name}'
    schema_pyproject_value = (
        f'nomad_simulation_parsers.schema_packages:{schema_entry_var_name}'
    )

    if dry_run:
        changed.extend(
            [
                str(parser_file),
                str(file_parser_file),
                str(parser_init),
                str(schema_file),
                str(parsers_registry),
                str(schema_registry),
                str(pyproject_path),
            ]
        )
        return changed

    parser_dir.mkdir(parents=True, exist_ok=True)
    parser_file.write_text(parser_content, encoding='utf-8')
    changed.append(str(parser_file))

    file_parser_file.write_text(file_parser_content, encoding='utf-8')
    changed.append(str(file_parser_file))

    if not parser_init.exists() or overwrite:
        parser_init.write_text('', encoding='utf-8')
        changed.append(str(parser_init))

    schema_file.write_text(schema_content, encoding='utf-8')
    changed.append(str(schema_file))

    anchor = f'{parser_entry_var_name} = EntryPoint('
    if _append_if_missing(parsers_registry, anchor, parser_entry_block):
        changed.append(str(parsers_registry))

    anchor = f'{schema_entry_var_name} = EntryPoint('
    if _append_if_missing(schema_registry, anchor, schema_entry_block):
        changed.append(str(schema_registry))

    added_entries = _insert_pyproject_entry_lines(
        pyproject_path,
        parser_entry_var_name,
        parser_pyproject_value,
        schema_entry_var_name,
        schema_pyproject_value,
    )
    if added_entries:
        changed.append(str(pyproject_path))

    return changed


def build_argument_parser(
    prog: str = 'nomad-sim-parser init',
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            'Initialize a parser and schema package from YAML metadata and register '
            'their plugin entry points.'
        ),
    )
    parser.add_argument(
        'metadata',
        help='Path to metadata YAML file.',
    )
    parser.add_argument(
        '--root',
        default='.',
        help='Path to the plugin repository root (default: current directory).',
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite generated parser/schema files if they already exist.',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print the files that would be touched without writing anything.',
    )
    parser.add_argument(
        '--use-mapping-parser',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Use mapping parsers in template (default: True).',
    )
    return parser


def run(
    arguments: list[str],
    *,
    prog: str = 'nomad-sim-parser init',
) -> int:
    """Parse initialization arguments and generate the requested files."""
    parser = build_argument_parser(prog=prog)
    args = parser.parse_args(arguments)

    metadata_path = Path(args.metadata).expanduser().resolve()
    root = Path(args.root).expanduser().resolve()

    if not metadata_path.exists():
        parser.error(f'Metadata file not found: {metadata_path}')

    try:
        changed_paths = initialize_from_metadata(
            metadata_path=metadata_path,
            root=root,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            use_mapping_parser=args.use_mapping_parser,
        )
    except Exception as exc:
        parser.exit(status=1, message=f'Error: {exc}\\n')

    mode_label = 'Would update' if args.dry_run else 'Updated'
    for path in changed_paths:
        print(f'{mode_label}: {path}')
    return 0
