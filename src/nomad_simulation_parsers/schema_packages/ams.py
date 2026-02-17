from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'ams_out'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.@')


class Representation(model_system.Representation):
    add_mapping_annotation(model_system.Representation.lattice_vectors, OUT_KEY, '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions, OUT_KEY, '.labels_positions[1]'
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def, OUT_KEY, '.labels_positions[0]'
    )
    add_mapping_annotation(
        model_system.Representation.m_def, OUT_KEY, '.lattice_vectors'
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, OUT_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        OUT_KEY,
        ('get_xc_functionals', ['.model_parameters.dft_potential']),
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.@')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(
        outputs.TotalEnergy.value, OUT_KEY, '.value || .energy_total'
    )
    add_mapping_annotation(
        outputs.TotalEnergy.contributions,
        OUT_KEY,
        ('get_energy_contributions', ['.energies']),
    )


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.value || .forces_total')
    add_mapping_annotation(
        outputs.TotalForce.contributions, OUT_KEY, ('get_contributions', ['.forces'])
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, OUT_KEY, '.eigenvalues')
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, OUT_KEY, '.occupations'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.@')
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        OUT_KEY,
        ('get_eigenvalues', ['.eigenvalues || .band_energy_ranges']),
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        model_method.DFT.m_def,
        OUT_KEY,
        '.geometry_optimization || molecular_dynamics || .single_point',
    )
    add_mapping_annotation(
        general.Simulation.model_system,
        OUT_KEY,
        '.geometry_optimization.step|| molecular_dynamics.step || .single_point',
    )
    add_mapping_annotation(
        general.Simulation.outputs,
        OUT_KEY,
        '.geometry_optimization.step|| molecular_dynamics.step || .single_point',
    )


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
