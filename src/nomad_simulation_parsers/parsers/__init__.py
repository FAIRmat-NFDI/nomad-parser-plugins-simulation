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
        description='Tolerance (in angstroms) for the cell positions to be considered '
        'equal.',
    )


abinit_parser = EntryPoint(
    name='parsers/abinit',
    aliases=['parsers/abinit'],
    description='NOMAD parser for ABINIT.',
    mainfile_contents_re=r'^\n*\.Version\s*[0-9.]*\s*of ABINIT\s*',
    python_package='nomad_simulation_parsers',
    parser_class_name='nomad_simulation_parsers.parsers.abinit.parser.AbinitParser',
    code_name='ABINIT',
    code_homepage='https://www.abinit.org/',
)

ams_parser = EntryPoint(
    name='parsers/ams',
    aliases=['parsers/ams'],
    description='NOMAD parser for AMS.',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'\* +\| +A M S +\| +\*',
    parser_class_name='nomad_simulation_parsers.parsers.ams.parser.AMSParser',
    code_name='AMS',
    code_homepage='https://www.scm.com',
)

crystal_parser = EntryPoint(
    name='parsers/crystal',
    aliases=['parsers/crystal'],
    description='NOMAD parser for CRYSTAL.',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'(\r?\n \*\s+CRYSTAL[\d]+\s+\*\r?\n \*\s*[a-zA-Z]+ :'
    r' \d+[\.\d+]*)',
    parser_class_name='nomad_simulation_parsers.parsers.crystal.parser.CrystalParser',
    code_name='CRYSTAL',
    code_homepage='https://www.crystal.unito.it/',
    code_category='Atomistic code',
)

exciting_parser = EntryPoint(
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

fhiaims_parser = EntryPoint(
    name='parsers/fhiaims',
    aliases=['parsers/fhi-aims', 'parsers/fhiaims'],
    description='NOMAD parser for FHIAIMS.',
    parser_class_name='nomad_simulation_parsers.parsers.fhiaims.parser.FHIAimsParser',
    python_package='nomad_simulation_parsers',
    code_name='fhiaims',
    code_homepage='https://aimsclub.fhi-berlin.mpg.de/',
    mainfile_contents_re=r'^(.*\n)*?\s*Invoking FHI-aims \.\.\.',
)

lammps_parser = EntryPoint(
    name='parsers/lammps',
    aliases=['parsers/lammps'],
    description='NOMAD parser for LAMMPS.',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'^LAMMPS\s+\(.+\)',
    parser_class_name='nomad_simulation_parsers.parsers.lammps.parser.LammpsParser',
    code_name='LAMMPS',
    code_homepage='https://lammps.sandia.gov/',
    code_category='Atomistic code',
)

gpaw_parser = EntryPoint(
    name='parsers/gpaw',
    aliases=['parsers/gpaw'],
    description='NOMAD parser for GPAW.',
    python_package='nomad_simulation_parsers',
    mainfile_mime_re='application/(x-tar|octet-stream)',
    mainfile_name_re=r'^.*\.(gpw2|gpw)$',
    parser_class_name='nomad_simulation_parsers.parsers.gpaw.parser.GPAWParser',
    code_name='GPAW',
    code_homepage='https://wiki.fysik.dtu.dk/gpaw/',
    code_category='Atomistic code',
)

h5md_parser = EntryPoint(
    name='parsers/h5md',
    aliases=['parsers/h5md'],
    description='NOMAD parser for H5MD.',
    python_package='nomad_simulation_parsers',
    mainfile_binary_header_re=b'^\\x89HDF',
    mainfile_contents_dict={'__has_all_keys': ['h5md']},
    mainfile_mime_re='(application/x-hdf)',
    mainfile_name_re=r'^.*\.(h5|hdf5)$',
    parser_class_name='nomad_simulation_parsers.parsers.h5md.parser.H5MDParser',
    code_name='H5MD',
    code_category='MD code',
    # metadata={
    #     'codeCategory': 'MD code',
    #     'codeLabel': 'H5MD',
    #     'codeLabelStyle': 'All in capitals',
    #     'codeName': 'h5md',
    #     # 'parserDirName': 'dependencies/parsers/atomistic/atomisticparsers/h5md/',
    #     'parserGitUrl': 'https://github.com/FAIRmat-NFDI/'
    #     'nomad-parser-plugins-simulation.git',  # ? Is this useful?
    #     # 'parserSpecific': '',
    #     # 'preamble': '',
    #     # 'status': 'production',
    #     # 'tableOfFiles': '',
    # },
)

octopus_parser = EntryPoint(
    name='parsers/octopus',
    aliases=['parsers/octopus'],
    description='NOMAD parser for OCTOPUS.',
    parser_class_name='nomad_simulation_parsers.parsers.octopus.parser.OctopusParser',
    python_package='nomad_simulation_parsers',
    mainfile_contents_re=r'\|0\) ~ \(0\) \|',
    code_name='Octopus',
    code_homepage='https://octopus-code.org/',
)

phonopy_parser = EntryPoint(
    name='parsers/phonopy',
    aliases=['parsers/phonopy'],
    description='NOMAD parser for PHONOPY.',
    mainfile_name_re='.*/phon[^/]+yaml',
    parser_class_name='nomad_simulation_parsers.parsers.phonopy.parser.PhonopyParser',
    code_name='phonopy',
    python_package='nomad_simulation_parsers',
    code_homepage='https://phonopy.github.io/phonopy/',
)

turbomole_parser = EntryPoint(
    name='parsers/turbomole',
    aliases=['parsers/turbomole'],
    description='NOMAD parser for TURBOMOLE.',
    mainfile_name_re='.*/Turbomole[^/]+',
    parser_class_name='nomad_simulation_parsers.parsers.turbomole.parser.TurbomoleParser',
    code_name='turbomole',
    python_package='nomad_simulation_parsers',
    code_homepage='https://www.turbomole.org/',
)

quantumespresso_parser = EntryPoint(
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

vasp_parser = EntryPoint(
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

wannier90_parser = Wannier90ParserEntryPoint(
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
