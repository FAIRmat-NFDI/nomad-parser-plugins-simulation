from nomad.config.models.plugins import ParserEntryPoint


class EntryPoint(ParserEntryPoint):
    def load(self):
        from nomad.parsing.parser import MatchingParserInterface

        return MatchingParserInterface(
            parser_class_name='nomad_simulation_parsers.parsers.exciting.parser.ExcitingParser',
            **self.dict(),
        )


exciting_parser_entry_point = EntryPoint(
    name='parsers/exciting',
    aliases=['parsers/exciting'],
    description='NOMAD parser for EXCITING.',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'EXCITING.*started[\s\S]+?All units are atomic ',
    mainfile_name_re=r'^.*.OUT(\.[^/]*)?$',
    code_name='exciting',
    code_homepage='http://exciting-code.org/',
)
