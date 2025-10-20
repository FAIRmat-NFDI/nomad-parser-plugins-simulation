from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations

m_package = SchemaPackage()

OUT_KEY = 'out'
INFO_KEY = 'info'
EIGENVALUES_KEY = 'eigenvalues'


class Program(general.Program):
    add_mapping_annotations(general.Program.version, OUT_KEY, '.Version')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotations(model_method.XCFunctional.libxc_name, OUT_KEY, '.@')


class DFT(model_method.DFT):
    add_mapping_annotations(
        model_method.DFT.xc_functionals,
        OUT_KEY,
        ('get_xc_functionals', ['.theory_level']),
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotations(model_system.ModelSystem.positions, OUT_KEY, '.positions')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotations(outputs.TotalEnergy.value, OUT_KEY, '.energy')
    add_mapping_annotations(outputs.TotalEnergy.value, INFO_KEY, '.Total || .value')


class TotalForce(outputs.TotalForce):
    add_mapping_annotations(outputs.TotalForce.value, INFO_KEY, '.forces')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotations(
        outputs.ElectronicEigenvalues.value, INFO_KEY, '.eigenvalues'
    )
    add_mapping_annotations(
        outputs.ElectronicEigenvalues.value, EIGENVALUES_KEY, '.eigenvalues'
    )
    add_mapping_annotations(
        outputs.ElectronicEigenvalues.occupation, INFO_KEY, '.occupations'
    )
    add_mapping_annotations(
        outputs.ElectronicEigenvalues.occupation, EIGENVALUES_KEY, '.occupations'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotations(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotations(outputs.Outputs.total_energies, INFO_KEY, '.energies')
    add_mapping_annotations(outputs.Outputs.total_forces, INFO_KEY, '.@')
    add_mapping_annotations(
        outputs.Outputs.electronic_eigenvalues,
        INFO_KEY,
        ('get_eigenvalues', ['eigenvalues.eigenvalues']),
    )
    add_mapping_annotations(
        outputs.Outputs.electronic_eigenvalues,
        EIGENVALUES_KEY,
        ('get_eigenvalues', ['eigenvalues']),
    )


class Simulation(general.Simulation):
    add_mapping_annotations(
        general.Simulation.program, OUT_KEY, ('get_header', ['.header'])
    )
    add_mapping_annotations(model_method.DFT.m_def, OUT_KEY, '.@')
    add_mapping_annotations(
        general.Simulation.model_system, OUT_KEY, ('get_systems', ['.minimization'])
    )
    add_mapping_annotations(
        general.Simulation.outputs,
        OUT_KEY,
        ('get_outputs', ['.minimization || .time_dependent[].iteration']),
    )
    add_mapping_annotations(general.Simulation.outputs, INFO_KEY, '.@')
    add_mapping_annotations(general.Simulation.outputs, EIGENVALUES_KEY, '.@')


add_mapping_annotations(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotations(general.Simulation.m_def, INFO_KEY, '@')
add_mapping_annotations(general.Simulation.m_def, EIGENVALUES_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
