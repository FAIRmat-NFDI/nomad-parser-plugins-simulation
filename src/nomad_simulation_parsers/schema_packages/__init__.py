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
