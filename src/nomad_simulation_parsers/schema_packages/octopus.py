from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'octopus_out'
INFO_KEY = 'octopus_info'
EIGENVALUES_KEY = 'octopus_eigenvalues'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.Version')


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, OUT_KEY, '.@')


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(out=Mapper(mapper='.@')))


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(out=Mapper(mapper=('get_xc_functionals', ['.theory_level']))))


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.value, OUT_KEY, '.energy')
    add_mapping_annotation(outputs.TotalEnergy.value, INFO_KEY, '.Total || .value')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, INFO_KEY, '.forces')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, INFO_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, EIGENVALUES_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, INFO_KEY, '.occupations'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, EIGENVALUES_KEY, '.occupations'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_energies, INFO_KEY, '.energies')
    add_mapping_annotation(outputs.Outputs.total_forces, INFO_KEY, '.@')
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        INFO_KEY,
        ('get_eigenvalues', ['eigenvalues.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        EIGENVALUES_KEY,
        ('get_eigenvalues', ['eigenvalues']),
    )


class Simulation(general.Simulation):
    add_mapping_annotation(
        general.Simulation.program, OUT_KEY, ('get_header', ['.header'])
    )
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.model_system, OUT_KEY, ('get_systems', ['.minimization'])
    )
    add_mapping_annotation(
        general.Simulation.outputs,
        OUT_KEY,
        ('get_outputs', ['.minimization || .time_dependent[].iteration']),
    )
    add_mapping_annotation(general.Simulation.outputs, INFO_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, EIGENVALUES_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, INFO_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EIGENVALUES_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
