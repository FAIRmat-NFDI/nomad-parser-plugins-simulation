from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field

from typing import Optional


class GSDParserEntryPoint(ParserEntryPoint):
    # parameter: int = Field(0, description='Custom configuration parameter')

    # def load(self):
    #     from nomad_parser_gsd.parsers.parser import GSDParser

    #     return GSDParser(**self.dict())
    parser_class_name: str = Field(
        description="""
        The fully qualified name of the Python class that implements the parser.
        This class must have a function `def parse(self, mainfile, archive, logger)`.
    """
    )
    code_name: Optional[str]
    code_homepage: Optional[str]
    code_category: Optional[str]
    metadata: Optional[dict] = Field(
        description="""
        Metadata passed to the UI. Deprecated. """
    )

    def load(self):
        from nomad.parsing import MatchingParserInterface

        return MatchingParserInterface(**self.dict())


parser_entry_point = GSDParserEntryPoint(
    name='nomad-parser-gsd',
    aliases=['gsd'],
    description='Parser for the GSD file format, the native file format for HOOMD-blue. (https://gsd.readthedocs.io/en/v3.3.1/, https://glotzerlab.engin.umich.edu/hoomd-blue/)',
    python_package='nomad-parser-gsd',
    mainfile_binary_header_re=b'^0x65DF65DF65DF65DF',
    # mainfile_contents_dict={},  # ?
    # mainfile_mime_re=,
    mainfile_name_re=r'^.*\.gsd$',
    parser_class_name='nomad_parser_gsd.parsers.parser.GSDParser',
    code_name='GSD',
    code_category='MD code',
    metadata={
        'codeCategory': 'MD code',
        'codeLabel': 'GSD',
        'codeLabelStyle': 'All in capitals',
        'codeName': 'gsd',
        'parserGitUrl': 'https://github.com/FAIRmat-NFDI/nomad-parser-gsd.git',
        'parserSpecific': '',
        'preamble': '',
        # 'status': 'production',
        'tableOfFiles': '',
    },
)
