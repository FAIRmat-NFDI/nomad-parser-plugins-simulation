from pathlib import Path

import pytest
import yaml

from nomad_simulation_parsers import mapping_report
from nomad_simulation_parsers.cli import main
from nomad_simulation_parsers.parser_init import initialize_from_metadata

PARSERS_INIT_TEMPLATE = """import importlib

from nomad.config.models.plugins import ParserEntryPoint
from nomad.utils import get_logger
from pydantic import Field

LOGGER = get_logger(__name__)


class EntryPoint(ParserEntryPoint):
    parser_class_name: str = Field(
        description='Parser class path.'
    )
"""


SCHEMA_INIT_TEMPLATE = """import importlib

from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class EntryPoint(SchemaPackageEntryPoint):
    module: str = Field(description='Module from which schema is loaded')
"""


PYPROJECT_TEMPLATE = """[project]
name = 'sample-plugin'

[project.entry-points.'nomad.plugin']
abinit_parser = 'nomad_simulation_parsers.parsers:abinit_parser'
abinit_schema_package = 'nomad_simulation_parsers.schema_packages:abinit_schema_package'

[tool.test]
"""


def _prepare_repo(tmp_path: Path) -> Path:
    root = tmp_path / 'repo'
    (root / 'src' / 'nomad_simulation_parsers' / 'parsers').mkdir(parents=True)
    (root / 'src' / 'nomad_simulation_parsers' / 'schema_packages').mkdir(parents=True)

    (root / 'src' / 'nomad_simulation_parsers' / 'parsers' / '__init__.py').write_text(
        PARSERS_INIT_TEMPLATE,
        encoding='utf-8',
    )
    (
        root / 'src' / 'nomad_simulation_parsers' / 'schema_packages' / '__init__.py'
    ).write_text(
        SCHEMA_INIT_TEMPLATE,
        encoding='utf-8',
    )
    (root / 'pyproject.toml').write_text(PYPROJECT_TEMPLATE, encoding='utf-8')
    return root


def test_unified_cli_dispatches_mapping_report(monkeypatch):
    called = {}

    def fake_report(arguments):
        called['arguments'] = arguments
        return 0

    monkeypatch.setattr(mapping_report, 'run', fake_report)

    assert main(['mapping-report', '--output', 'report.md']) == 0
    assert called['arguments'] == ['--output', 'report.md']


def test_unified_cli_visualizes_parsed_blocks(tmp_path: Path, capsys, monkeypatch):
    mainfile = Path(__file__).parent / 'data' / 'exciting' / 'C_minimal' / 'INFO.OUT'
    output = tmp_path / 'parsed-blocks.html'
    opened = []
    monkeypatch.setattr('webbrowser.open_new_tab', opened.append)

    assert (
        main(
            [
                'view-blocks',
                'InfoFileParser',
                str(mainfile),
                '--output',
                str(output),
            ]
        )
        == 0
    )

    assert capsys.readouterr().out.strip() == str(output.resolve())
    assert opened == [output.resolve().as_uri()]
    html = output.read_text(encoding='utf-8')
    assert '<mark' in html


def test_unified_cli_requires_a_subcommand():
    assert main(['--help']) == 0


def test_unified_cli_rejects_unknown_command():
    with pytest.raises(SystemExit) as error:
        main(['unknown-command'])

    assert error.value.code == 2


def test_initialize_from_metadata_creates_parser_and_schema(tmp_path: Path):
    root = _prepare_repo(tmp_path)

    metadata = {
        'code_key': 'my-code',
        'name': 'MyCode',
        'description': 'NOMAD parser for MyCode.',
        'code_homepage': 'https://example.org/mycode',
        'parser': {
            'mainfile_name_re': r'^.*\\.mycode$',
            'aliases': ['parsers/my-code', 'parsers/mycode'],
        },
    }
    metadata_path = root / 'parser_metadata.yaml'
    metadata_path.write_text(yaml.safe_dump(metadata), encoding='utf-8')

    changed_paths = initialize_from_metadata(metadata_path=metadata_path, root=root)

    parser_file = (
        root / 'src' / 'nomad_simulation_parsers' / 'parsers' / 'my_code' / 'parser.py'
    )
    file_parser_file = (
        root
        / 'src'
        / 'nomad_simulation_parsers'
        / 'parsers'
        / 'my_code'
        / 'file_parser.py'
    )
    schema_file = (
        root / 'src' / 'nomad_simulation_parsers' / 'schema_packages' / 'my_code.py'
    )

    assert parser_file.exists()
    assert file_parser_file.exists()
    assert schema_file.exists()

    parser_content = parser_file.read_text(encoding='utf-8')
    assert 'class MyCodeParser(MatchingParser):' in parser_content
    assert 'class MyCodeParserArchiveWriter(ArchiveWriter):' in parser_content
    assert 'class MyCodeParserMainfileParser(TextParser):' in parser_content
    assert 'class MyCodeParserMetainfoParser(MetainfoParser):' in parser_content
    assert 'Simulation(program=Program(name=self.code_name))' in parser_content
    assert 'from .file_parser import MyCodeParserOutParser' in parser_content

    file_parser_content = file_parser_file.read_text(encoding='utf-8')
    assert 'class MyCodeParserOutParser(TextParser):' in file_parser_content

    schema_content = schema_file.read_text(encoding='utf-8')
    assert "OUT_KEY = 'my_code'" in schema_content
    assert (
        'add_mapping_annotation(general.Simulation.m_def, OUT_KEY, ' in schema_content
    )

    parsers_init = (
        root / 'src' / 'nomad_simulation_parsers' / 'parsers' / '__init__.py'
    ).read_text(encoding='utf-8')
    assert 'my_code_parser = EntryPoint(' in parsers_init
    assert 'mainfile_name_re=' in parsers_init
    assert 'mycode' in parsers_init

    schema_init = (
        root / 'src' / 'nomad_simulation_parsers' / 'schema_packages' / '__init__.py'
    ).read_text(encoding='utf-8')
    assert 'my_code_schema_package = EntryPoint(' in schema_init

    pyproject = (root / 'pyproject.toml').read_text(encoding='utf-8')
    assert (
        'my_code_parser = "nomad_simulation_parsers.parsers:my_code_parser"'
        in pyproject
    )
    assert (
        'my_code_schema_package = "'
        'nomad_simulation_parsers.schema_packages:my_code_schema_package"' in pyproject
    )

    assert str(parser_file) in changed_paths
    assert str(file_parser_file) in changed_paths
    assert str(schema_file) in changed_paths


def test_initialize_from_metadata_without_mapping_parser(tmp_path: Path):
    root = _prepare_repo(tmp_path)

    metadata = {
        'code_key': 'plain-parser',
        'name': 'PlainParser',
        'description': 'NOMAD parser for PlainParser.',
        'parser': {
            'use_mapping_parser': False,
        },
    }
    metadata_path = root / 'parser_metadata.yaml'
    metadata_path.write_text(yaml.safe_dump(metadata), encoding='utf-8')

    initialize_from_metadata(
        metadata_path=metadata_path, root=root, use_mapping_parser=False
    )

    parser_file = (
        root
        / 'src'
        / 'nomad_simulation_parsers'
        / 'parsers'
        / 'plain_parser'
        / 'parser.py'
    )
    schema_file = (
        root
        / 'src'
        / 'nomad_simulation_parsers'
        / 'schema_packages'
        / 'plain_parser.py'
    )

    parser_content = parser_file.read_text(encoding='utf-8')
    assert 'class PlainParserParserArchiveWriter(ArchiveWriter):' in parser_content
    assert 'class PlainParserParserMainfileParser(TextParser):' not in parser_content
    assert (
        'class PlainParserParserMetainfoParser(MetainfoParser):' not in parser_content
    )
    assert 'mainfile_parser = PlainParserParserOutParser()' in parser_content
    assert 'self.mainfile_parser.parse()' in parser_content

    schema_content = schema_file.read_text(encoding='utf-8')
    assert 'add_mapping_annotation' not in schema_content
    assert 'OUT_KEY =' not in schema_content


def test_initialize_from_metadata_dry_run_does_not_write_files(tmp_path: Path):
    root = _prepare_repo(tmp_path)

    metadata = {
        'code_key': 'dryrun',
        'name': 'DryRun',
        'description': 'NOMAD parser for DryRun.',
    }
    metadata_path = root / 'parser_metadata.yaml'
    metadata_path.write_text(yaml.safe_dump(metadata), encoding='utf-8')

    changed_paths = initialize_from_metadata(
        metadata_path=metadata_path,
        root=root,
        dry_run=True,
    )

    parser_file = (
        root / 'src' / 'nomad_simulation_parsers' / 'parsers' / 'dryrun' / 'parser.py'
    )
    schema_file = (
        root / 'src' / 'nomad_simulation_parsers' / 'schema_packages' / 'dryrun.py'
    )

    file_parser_file = (
        root
        / 'src'
        / 'nomad_simulation_parsers'
        / 'parsers'
        / 'dryrun'
        / 'file_parser.py'
    )

    assert not parser_file.exists()
    assert not file_parser_file.exists()
    assert not schema_file.exists()
    assert any(path.endswith('pyproject.toml') for path in changed_paths)
