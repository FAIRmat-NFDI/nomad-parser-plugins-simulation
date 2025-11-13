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
        schema_package.MagneticShielding.value, 
        GIPAW_XML_KEY, 
        ('get_magnetic_shieldings', ['.@'])
    )


class MagneticSusceptibility(schema_package.MagneticSusceptibility):
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value, 
        GIPAW_XML_KEY, 
        ('get_magnetic_susceptibilities', ['.@'], dict(name='value'))
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value_vgv_approx, 
        GIPAW_OUT_KEY, 
        '.value_vgv_approx'
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value_vgv_approx, 
        GIPAW_XML_KEY, 
        ('get_magnetic_susceptibilities', ['.@'], dict(name='susceptibility_low'))
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value_pgv_approx, 
        GIPAW_OUT_KEY, 
        '.value_pgv_approx'
    )
    add_mapping_annotation(
        schema_package.MagneticSusceptibility.value_pgv_approx, 
        GIPAW_XML_KEY, 
        ('get_magnetic_susceptibilities', ['.@'], dict(name='susceptibility_high'))
    )


class ElectricFieldGradient(schema_package.ElectricFieldGradient):
    add_mapping_annotation(
        schema_package.ElectricFieldGradient.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.ElectricFieldGradient.value, GIPAW_XML_KEY, ('get_efg', ['.@'])
    )


class HyperfineDipolar(schema_package.HyperfineDipolar):
    add_mapping_annotation(
        schema_package.HyperfineDipolar.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.HyperfineDipolar.value, 
        GIPAW_XML_KEY, 
        ('get_hyperfine_dipolar', ['.@'])
    )


class HyperfineFermiContact(schema_package.HyperfineFermiContact):
    add_mapping_annotation(
        schema_package.HyperfineFermiContact.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.HyperfineFermiContact.value, 
        GIPAW_XML_KEY, 
        ('get_hyperfine_fermi_contact', ['.@'])
    )


class DeltaG(schema_package.DeltaG):
    add_mapping_annotation(
        schema_package.DeltaG.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.DeltaG.value, GIPAW_XML_KEY, ('get_delta_g', ['.@'])
    )


class DeltaGParatec(schema_package.DeltaGParatec):
    add_mapping_annotation(
        schema_package.DeltaGParatec.value, GIPAW_OUT_KEY, '.value'
    )
    add_mapping_annotation(
        schema_package.DeltaGParatec.value, 
        GIPAW_XML_KEY, 
        ('get_delta_g', ['.@'])
    )


class Outputs(schema_package.Outputs):
    add_mapping_annotation(
        schema_package.Outputs.m_def, GIPAW_OUT_KEY, ('get_gipaw_text', ['.@']),
    )
    add_mapping_annotation(schema_package.Outputs.m_def, GIPAW_XML_KEY, '.output')
    add_mapping_annotation(
        schema_package.Outputs.magnetic_shieldings, 
        GIPAW_OUT_KEY, 
        '.magnetic_shieldings'
    )
    add_mapping_annotation(
        schema_package.Outputs.magnetic_shieldings, 
        GIPAW_XML_KEY, 
        '.shielding_tensors.atom'
    )
    add_mapping_annotation(
        schema_package.Outputs.magnetic_susceptibilities, 
        GIPAW_OUT_KEY, 
        '.magnetic_susceptibilities'
    )
    add_mapping_annotation(
        schema_package.Outputs.magnetic_susceptibilities, GIPAW_XML_KEY, '.@'
    )
    add_mapping_annotation(
        schema_package.Outputs.electric_field_gradients, 
        GIPAW_OUT_KEY, 
        '.electric_field_gradients'
    )
    add_mapping_annotation(
        schema_package.Outputs.electric_field_gradients, 
        GIPAW_XML_KEY, 
        '.electric_field_gradients.atom'
    )
    add_mapping_annotation(
        schema_package.Outputs.hyperfine_dipolar, GIPAW_OUT_KEY, '.hyperfine_dipolar'
    )
    add_mapping_annotation(
        schema_package.Outputs.hyperfine_dipolar, 
        GIPAW_XML_KEY, 
        '.hyperfine_dipolar.atom'
    )
    add_mapping_annotation(
        schema_package.Outputs.hyperfine_fermi_contact, 
        GIPAW_OUT_KEY, 
        '.hyperfine_fermi_contact'
    )
    add_mapping_annotation(
        schema_package.Outputs.hyperfine_fermi_contact, 
        GIPAW_XML_KEY, 
        '.hyperfine_fermi_contact.atom'
    )
    add_mapping_annotation(
        schema_package.Outputs.delta_g_paratec, GIPAW_OUT_KEY, '.delta_g_paratec'
    )
    add_mapping_annotation(
        schema_package.Outputs.delta_g_paratec, GIPAW_XML_KEY, '.delta_g_paratec'
    )
    add_mapping_annotation(
        schema_package.Outputs.delta_g, GIPAW_OUT_KEY, '.delta_g'
    )
    add_mapping_annotation(
        schema_package.Outputs.delta_g, GIPAW_XML_KEY, '.delta_g'
    )


add_mapping_annotation(general.Simulation.m_def, GIPAW_OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, GIPAW_XML_KEY, '@')

class Simulation(general.Simulation):
    pass


try:
    m_package.__init_metainfo__()
except Exception:
    pass
