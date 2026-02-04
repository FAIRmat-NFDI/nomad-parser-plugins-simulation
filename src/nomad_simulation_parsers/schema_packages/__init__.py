import importlib

from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class EntryPoint(SchemaPackageEntryPoint):
    module: str = Field(description="""Module from which schema is loaded""")

    def load(self):
        try:
            return importlib.import_module(self.module).m_package
        except Exception:
            return None


# in schema_packages/__init__.py
# from nomad.utils import get_logger

# LOGGER = get_logger(__name__)


# class EntryPoint(SchemaPackageEntryPoint):
#     module: str = Field(description='Module from which schema is loaded')

#     def load(self):
#         try:
#             mod = importlib.import_module(self.module)
#             mp = getattr(mod, 'm_package', None)
#             if mp is None:
#                 raise AttributeError(f'{self.module} has no m_package')
#             return mp
#         except Exception as e:
#             LOGGER.error(f'Could not load schema package {self.module}', exc_info=e)
#             return None


abinit_schema_package = EntryPoint(
    name='AbinitSchemaPackage',
    description='Schema package for abinit.',
    module='nomad_simulation_parsers.schema_packages.abinit',
)


ams_schema_package = EntryPoint(
    name='AMSSchemaPackage',
    description='Schema package for AMS.',
    module='nomad_simulation_parsers.schema_packages.ams',
)

crystal_schema_package = EntryPoint(
    name='CrystalSchemaPackage',
    description='Schema package for Crystal.',
    module='nomad_simulation_parsers.schema_packages.crystal',
)

exciting_schema_package = EntryPoint(
    name='ExcitingSchemaPackage',
    description='Schema package for exciting.',
    module='nomad_simulation_parsers.schema_packages.exciting',
)

fhiaims_schema_package = EntryPoint(
    name='FHIAimsSchemaPackage',
    description='Schema package for FHIAims.',
    module='nomad_simulation_parsers.schema_packages.fhiaims',
)

gaussian_schema_package = EntryPoint(
    name='GaussianSchemaPackage',
    description='Schema package for Gaussian.',
    module='nomad_simulation_parsers.schema_packages.gaussian',
)

gpaw_schema_package = EntryPoint(
    name='GPAWSchemaPackage',
    description='Schema package for GPAW.',
    module='nomad_simulation_parsers.schema_packages.gpaw',
)

gromacs_schema_package = EntryPoint(
    name='GromacsSchemaPackage',
    description='Schema package for Gromacs.',
    module='nomad_simulation_parsers.schema_packages.gromacs',
)

h5md_schema_package = EntryPoint(
    name='H5MDSchemaPackage',
    description='Schema package for H5MD.',
    module='nomad_simulation_parsers.schema_packages.h5md',
)

octopus_schema_package = EntryPoint(
    name='OctopusSchemaPackage',
    description='Schema package for Octopus.',
    module='nomad_simulation_parsers.schema_packages.octopus',
)

phonopy_schema_package = EntryPoint(
    name='PhonopySchemaPackage',
    description='Schema package for Phonopy.',
    module='nomad_simulation_parsers.schema_packages.phonopy',
)

quantumespresso_schema_package = EntryPoint(
    name='QuantumEspressoSchemaPackage',
    description='Schema package for Quantum Espresso.',
    module='nomad_simulation_parsers.schema_packages.quantumespresso.common',
)

vasp_schema_package = EntryPoint(
    name='VASPSchemaPackage',
    description='Schema package for VASP.',
    module='nomad_simulation_parsers.schema_packages.vasp',
)

wannier90_schema_package = EntryPoint(
    name='Wannier90SchemaPackage',
    description='Schema package for Wannier90.',
    module='nomad_simulation_parsers.schema_packages.wannier90',
)

yambo_schema_package = EntryPoint(
    name='YamboSchemaPackage',
    description='Schema package for Yambo.',
    module='nomad_simulation_parsers.schema_packages.yambo',
)
