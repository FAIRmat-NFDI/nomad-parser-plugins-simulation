import importlib

from nomad.config.models.plugins import ParserEntryPoint
from nomad.utils import get_logger
from pydantic import Field

LOGGER = get_logger(__name__)


class EntryPoint(ParserEntryPoint):
    parser_class_name: str = Field(
        description="""
        The fully qualified name of the Python class that implements the parser.
        This class must have a function `def parse(self, mainfile, archive, logger)`.
    """
    )

    def load(self):
        try:
            module_path, cls_name = self.parser_class_name.rsplit('.', 1)
            module = importlib.import_module(module_path)
            cls = getattr(module, cls_name)
            return cls(**self.model_dump(exclude={'parser_class_name'}))
        except Exception as e:
            LOGGER.error(
                f'Could not load parser class {self.parser_class_name}', exc_info=e
            )


class Wannier90ParserEntryPoint(EntryPoint):
    equal_cell_positions_tolerance: float = Field(
        1e-2,
        description='Tolerance (in angstroms) for the cell positions to be considered'
        'equal.',
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

fhiaims_parser_entry_point = EntryPoint(
    name='parsers/fhiaims',
    aliases=['parsers/fhi-aims', 'parsers/fhiaims'],
    description='NOMAD parser for FHIAIMS.',
    parser_class_name='nomad_simulation_parsers.parsers.fhiaims.parser.FHIAimsParser',
    python_package='nomad_simulation_parsers',
    code_name='phonopy',
    code_homepage='https://aimsclub.fhi-berlin.mpg.de/',
    mainfile_contents_re=r'^(.*\n)*?\s*Invoking FHI-aims \.\.\.',
)

phonopy_parser_entry_point = EntryPoint(
    name='parsers/phonopy',
    aliases=['parsers/phonopy'],
    description='NOMAD parser for PHONOPY.',
    mainfile_name_re='.*/phon[^/]+yaml',
    parser_class_name='nomad_simulation_parsers.parsers.phonopy.parser.PhonopyParser',
    code_name='phonopy',
    python_package='nomad_simulation_parsers',
    code_homepage='https://phonopy.github.io/phonopy/',
)

quantumespresso_parser_entry_point = EntryPoint(
    name='parsers/quantumespresso',
    aliases=['parsers/quantumespresso'],
    description='NOMAD parser for QUANTUMESPRESSO.',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=(
        r'(Program [A-Z]+.*starts)|(Current dimensions of program [A-Z]+ are)'
    ),
    supported_compressions=['gz', 'bz2', 'xz'],
    parser_class_name='nomad_simulation_parsers.parsers.quantumespresso.parser.QuantumEspressoParser',
    code_name='QuantumESPRESSO',
    code_homepage='http://www.quantum-espresso.org/',
    code_category='Atomistic code',
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

wannier90_parser_entry_point = Wannier90ParserEntryPoint(
    name='parsers/wannier90',
    aliases=['parsers/wannier90'],
    description='NOMAD parser for WANNIER90.',
    parser_class_name='nomad_simulation_parsers.parsers.wannier90.parser.Wannier90Parser',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'\|\s*WANNIER90\s*\|',
    code_name='Wannier90',
    code_homepage='http://www.wannier.org/',
    code_category='Atomistic code',
)
