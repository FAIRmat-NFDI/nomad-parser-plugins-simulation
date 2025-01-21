from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class EntryPoint(ParserEntryPoint):
    parser_class_name: str = Field(
        description="""
        The fully qualified name of the Python class that implements the parser.
        This class must have a function `def parse(self, mainfile, archive, logger)`.
    """
    )

    def load(self):
        from nomad.parsing.parser import MatchingParserInterface

        return MatchingParserInterface(
            self.parser_class_name,
            **self.dict(),
        )


exciting_parser_entry_point = EntryPoint(
    name='parsers/exciting',
    aliases=['parsers/exciting'],
    description='NOMAD parser for EXCITING.',
    parser_class_name='nomad_simulation_parsers.parsers.exciting.parser.ExcitingParser',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'EXCITING.*started[\s\S]+?All units are atomic ',
    mainfile_name_re=r'^.*.OUT(\.[^/]*)?$',
    code_name='exciting',
    code_homepage='http://exciting-code.org/',
)

vasp_parser_entry_point = EntryPoint(
    name='parsers/vasp',
    description='Parser for VASP XML and OUTCAR outputs',
    parser_class_name='nomad_simulation_parsers.parsers.vasp.parser.VASPParser',
    python_package='nomad_simulation_parsers',
    code_name='VASP',
    mainfile_contents_re=(
        r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*?\s*<modeling>?\s*'
        r'<generator>?\s*<i name="program" type="string">\s*vasp\s*</i>?|'
        r'^\svasp[\.\d]+.+?(?:\(build|complex)[\s\S]+?executed on'
    ),
    mainfile_mime_re='(application/.*)|(text/.*)',
    mainfile_name_re='.*[^/]*xml[^/]*',
    mainfile_alternative=True,
    supported_compressions=['gz', 'bz2', 'xz'],
)
