from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general, model_system, outputs

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()
OUT_KEY = 'cp2k_out'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, ('get_version', ['.@']))


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.label')
    add_mapping_annotation(model_system.AtomsState.atomic_number, OUT_KEY, '.number')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, OUT_KEY, '.lattice_vectors'
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, OUT_KEY, '.@')
    add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, OUT_KEY, '.atoms')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.value, OUT_KEY, '.energy')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.forces')


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(outputs.SCFSteps.energies_total, OUT_KEY, '.energies_total')
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, OUT_KEY, '.delta_energies_total'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.scf_steps, OUT_KEY, '.scf_steps')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.model_system, OUT_KEY, ('get_systems', ['.@'])
    )
    add_mapping_annotation(general.Simulation.outputs, OUT_KEY, ('get_outputs', ['.@']))


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
m_package.__init_metainfo__()
