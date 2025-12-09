import json
import os
from collections.abc import Generator
from contextlib import contextmanager

import pytest
from nomad.datamodel import EntryArchive

# from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers import abinit_parser
from nomad_simulation_parsers.parsers.abinit.file_parser import AbinitOutParser

from ..conftest import compare_values, get_child_archive_keys

LOGGER = get_logger(__name__)
PARSERS = {}


def join_test_dir(*path_segments) -> str:
    return os.path.join('tests/data/abinit', *path_segments)


@pytest.fixture(scope='session')
def parser() -> MatchingParser:
    return abinit_parser.load()


@contextmanager
def generate_out_parser(mainfile: str) -> Generator[AbinitOutParser, None, None]:
    parser = PARSERS.get(mainfile, {}).get('out')
    if parser is None:
        parser = AbinitOutParser()
        parser.mainfile = mainfile
        parser.parse()
        PARSERS.setdefault(mainfile, {})['out'] = parser
    try:
        yield parser
    finally:
        parser.close()


class TestAbinitParser:
    @pytest.mark.parametrize(
        'filename, expected_keys',
        [('Fe/Fe.out', True), ('ZrO2_GW/A1.abo', ['GW', 'GW_workflow'])],
    )
    def test_is_mainfile(self, filename, expected_keys, parser):
        assert get_child_archive_keys(join_test_dir(filename), parser) == expected_keys


class TestAbinitOutParser:
    @pytest.mark.parametrize(
        'filename', [('Fe/Fe.out'), ('H2/H2.out'), ('ZrO2_GW/A1.abo')]
    )
    def test_quantities(self, filename):
        with (
            generate_out_parser(join_test_dir(filename)) as file_parser,
            open(
                join_test_dir(
                    os.path.dirname(filename),
                    f'reference_{os.path.basename(filename)}.json',
                )
            ) as reference_file,
        ):
            assert compare_values(
                file_parser.results, json.load(reference_file), all_keys=True
            )

    @pytest.mark.parametrize(
        'filename, reference_input_vars', [('Fe/Fe.out', {'znucl': [26.0, 26.0]})]
    )
    def test_input_vars(self, filename, reference_input_vars):
        with generate_out_parser(join_test_dir(filename)) as file_parser:
            assert compare_values(file_parser.input_vars, reference_input_vars)


def test_parse_file(parser):
    archive = EntryArchive()
    parser.parse('tests/data/abinit/Fe/Fe.out', archive, LOGGER)
