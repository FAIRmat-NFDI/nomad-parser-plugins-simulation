from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
    workflow,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'out'
DOS_KEY = 'dos'


class GeometryOptimizationMethod(
    workflow.geometry_optimization.GeometryOptimizationMethod
):
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.optimization_method,
        OUT_KEY,
        ('get_workflow_method', []),
    )
"""    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.convergence_tolerance_energy_difference,
        OUT_KEY,
        ('get_input_var', [], dict(name='tolmxde', n_dataset=1, default=0.0)),
        unit='hartree',
    )
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.convergence_tolerance_force_maximum,
        OUT_KEY,
        ('get_input_var', [], dict(name='tolmxf', n_dataset=1, default=0.0)),
        unit='hartree/bohr',
    )
"""
# Geometry Optimization

add_mapping_annotation(
    workflow.geometry_optimization.GeometryOptimizationMethod.m_def, OUT_KEY, '@'
)
add_mapping_annotation(workflow.GeometryOptimization.m_def, OUT_KEY, '@')

class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.label')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors,
        OUT_KEY,
        'dataset[0].x_abinit_vprim',
        unit='bohr',
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, OUT_KEY, '.@')
    add_mapping_annotation(
        model_system.AtomsState.m_def, OUT_KEY, ('get_atoms', []), cache=True
    )
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        OUT_KEY,
        '.cartesian_coordinates',
        unit='bohr',
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(
        model_method.XCComponent.canonical_label, OUT_KEY, '.XC_functional_name'
    )


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components, OUT_KEY, ('get_xc_functionals', [])
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.@')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(
        outputs.TotalEnergy.value, OUT_KEY, '.value || .energy_total'
    )
    add_mapping_annotation(outputs.TotalEnergy.name, OUT_KEY, '.name')
    add_mapping_annotation(
        outputs.TotalEnergy.contributions, OUT_KEY, ('get_energy_contributions', ['.@'])
    )


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.cartesian_forces')


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value, OUT_KEY, '.value', unit='1 / hartree'
    )


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    add_mapping_annotation(
        outputs.ElectronicBandStructure.value, OUT_KEY, '.energies', unit='hartree'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.@')
    add_mapping_annotation(
        outputs.Outputs.electronic_dos, DOS_KEY, ('get_dos', ['.data'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        OUT_KEY,
        ('get_bandstructures', ['.eigenvalues', '.occupation_numbers']),
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.datetime,
        OUT_KEY,
        ('get_datetime', ['x_abinit_start_date', 'x_abinit_start_time']),
    )
    add_mapping_annotation(
        general.Simulation.model_system, OUT_KEY, ('get_systems', [])
    )
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, OUT_KEY, ('get_outputs', []))
    add_mapping_annotation(general.Simulation.outputs, DOS_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOS_KEY, '.@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
