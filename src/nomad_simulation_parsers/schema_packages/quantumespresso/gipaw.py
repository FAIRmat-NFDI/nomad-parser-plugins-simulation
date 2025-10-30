from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations, remove_mapping_annotations

from .common import OUT_KEY, XML_KEY

from nomad_nmr_schema.schema_packages import schema_package

m_package = SchemaPackage()


class MagneticShielding(schema_package.MagneticShielding):
    add_mapping_annotations(
        schema_package.MagneticShielding.value, OUT_KEY, '.@'
    )


class Outputs(schema_package.Outputs):
    remove_mapping_annotations(general.Simulation.outputs)
    add_mapping_annotations(
        schema_package.Outputs.m_def,
        OUT_KEY,
        ('get_magnetic_shieldings', ['.@']),
    )
    add_mapping_annotations(schema_package.Outputs.magnetic_shieldings, OUT_KEY, '.@')


class Simulation(general.Simulation):
    pass


try:
    m_package.__init_metainfo__()
except Exception:
    pass
