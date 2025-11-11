from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

from .common import OUT_KEY, XML_KEY, GIPAW_OUT_KEY, GIPAW_XML_KEY

from nomad_nmr_schema.schema_packages import schema_package

m_package = SchemaPackage()


class MagneticShielding(schema_package.MagneticShielding):
    add_mapping_annotation(
        schema_package.MagneticShielding.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.MagneticShielding.value, GIPAW_XML_KEY, ('get_nmr_xml', ['.@'], dict(name='pippo'))
    )


class MagneticSusceptibility(schema_package.MagneticSusceptibility):

    # txt = ''
    # for item in dir(schema_package.MagneticSusceptibility):
    #     txt += f'{item}\n'
    # with open('SUS.txt', 'w') as f:
    #     f.write(txt)

    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value_vgv_approx, GIPAW_OUT_KEY, '.value_vgv_approx'
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value_pgv_approx, GIPAW_OUT_KEY, '.value_pgv_approx'
    )


class ElectricFieldGradient(schema_package.ElectricFieldGradient):

    # txt = ''
    # for item in dir(schema_package.ElectricFieldGradient):
    #     txt += f'{item}\n'
    # with open('EFG.txt', 'w') as f:
    #     f.write(txt)

    add_mapping_annotation(
        schema_package.ElectricFieldGradient.value, GIPAW_OUT_KEY, '.value'
    )


class Outputs(schema_package.Outputs):
    add_mapping_annotation(
        schema_package.Outputs.m_def,
        GIPAW_OUT_KEY,
        ('get_nmr_text', ['.@']),
    )
    add_mapping_annotation(
        schema_package.Outputs.m_def,
        GIPAW_XML_KEY,
        '@',
    )
    add_mapping_annotation(schema_package.Outputs.magnetic_shieldings, GIPAW_OUT_KEY, '.magnetic_shieldings')
    add_mapping_annotation(schema_package.Outputs.magnetic_shieldings, GIPAW_XML_KEY, '.magnetic_shieldings')
    add_mapping_annotation(schema_package.Outputs.magnetic_susceptibilities, GIPAW_OUT_KEY, '.magnetic_susceptibilities')
    add_mapping_annotation(schema_package.Outputs.electric_field_gradients, GIPAW_OUT_KEY, '.electric_field_gradients')

add_mapping_annotation(general.Simulation.m_def, GIPAW_OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, GIPAW_XML_KEY, '@')

class Simulation(general.Simulation):
    pass


try:
    m_package.__init_metainfo__()
except Exception:
    pass
