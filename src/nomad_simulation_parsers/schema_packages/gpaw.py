from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

GPW_KEY = 'gpw'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, GPW_KEY, '.program_version')


class AtomicCell(model_system.AtomicCell):
    add_mapping_annotation(
        model_system.AtomicCell.lattice_vectors, GPW_KEY, '.unitcell'
    )
    add_mapping_annotation(
        model_system.AtomicCell.periodic_boundary_conditions,
        GPW_KEY,
        '.boundary_conditions',
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, GPW_KEY, '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions, GPW_KEY, '.atom_positions'
    )
    add_mapping_annotation(model_system.AtomicCell.m_def, GPW_KEY, '.@')
    add_mapping_annotation(model_system.AtomsState.m_def, GPW_KEY, '.labels')


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, GPW_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components, GPW_KEY, '.xcfunctional'
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc.m_annotations, GPW_KEY, '.@')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.value, GPW_KEY, '.total || .value')
    add_mapping_annotation(outputs.TotalEnergy.name, GPW_KEY, '.name')
    add_mapping_annotation(outputs.TotalEnergy.contributions, GPW_KEY, '.contributions')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, GPW_KEY, '.total || .value')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, GPW_KEY, '.eigenvalues')
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, GPW_KEY, '.occupations'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(
        outputs.Outputs.total_energies, GPW_KEY, ('get_energies', [])
    )
    add_mapping_annotation(outputs.Outputs.total_forces, GPW_KEY, ('get_forces', []))
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues, GPW_KEY, ('get_eigenvalues', [])
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, GPW_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_system, GPW_KEY, '.@')
    add_mapping_annotation(model_method.DFT.m_def, GPW_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, GPW_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, GPW_KEY, '.@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
