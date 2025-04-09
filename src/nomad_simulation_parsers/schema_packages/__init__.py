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


exciting_schema_package_entry_point = EntryPoint(
    name='ExcitingSchemaPackage',
    description='Schema package for exciting.',
    module='nomad_simulation_parsers.schema_packages.exciting',
)

fhiaims_schema_package_entry_point = EntryPoint(
    name='FHIAimsSchemaPackage',
    description='Schema package for FHIAims.',
    module='nomad_simulation_parsers.schema_packages.fhiaims',
)

phonopy_schema_package_entry_point = EntryPoint(
    name='PhonopySchemaPackage',
    description='Schema package for Phonopy.',
    module='nomad_simulation_parsers.schema_packages.phonopy',
)

vasp_schema_package_entry_point = EntryPoint(
    name='VASPSchemaPackage',
    description='Schema package for VASP.',
    module='nomad_simulation_parsers.schema_packages.vasp',
)

wannier90_schema_package_entry_point = EntryPoint(
    name='Wannier90SchemaPackage',
    description='Schema package for Wannier90.',
    module='nomad_simulation_parsers.schema_packages.wannier90',
)
