from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations, remove_mapping_annotations

from .common import OUT_KEY, XML_KEY, GPAW_OUT_KEY

from nomad_nmr_schema.schema_packages import schema_package

m_package = SchemaPackage()


class MagneticShielding(schema_package.MagneticShielding):
    add_mapping_annotations(
        schema_package.MagneticShielding.value, GPAW_OUT_KEY, '.value'
    )


class MagneticSusceptibility(schema_package.MagneticSusceptibility):
    add_mapping_annotations(
        schema_package.MagneticSusceptibility.value, GPAW_OUT_KEY, '.value'
    )
    add_mapping_annotations(
        schema_package.MagneticSusceptibility.value_vgv_approx, GPAW_OUT_KEY, '.value_vgv_approx'
    )
    add_mapping_annotations(
        schema_package.MagneticSusceptibility.value_pgv_approx, GPAW_OUT_KEY, '.value_pgv_approx'
    )


class Outputs(schema_package.Outputs):
    add_mapping_annotations(
        schema_package.Outputs.m_def,
        GPAW_OUT_KEY,
        ('get_nmr', ['.@']),
    )
    add_mapping_annotations(schema_package.Outputs.magnetic_shieldings, GPAW_OUT_KEY, '.magnetic_shieldings')
    add_mapping_annotations(schema_package.Outputs.magnetic_susceptibilities, GPAW_OUT_KEY, '.magnetic_susceptibilities')

add_mapping_annotations(general.Simulation.m_def, GPAW_OUT_KEY, '@')

class Simulation(general.Simulation):
    pass


try:
    m_package.__init_metainfo__()
except Exception:
    pass
