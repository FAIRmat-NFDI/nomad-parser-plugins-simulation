from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

OUT_KEY = 'gaussian_out'

m_package = SchemaPackage()


class Program(general.Program):
    add_mapping_annotation(general.Program.name, OUT_KEY, '.name')
    add_mapping_annotation(general.Program.version, OUT_KEY, '.version')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        OUT_KEY,
        '.positions',
    )
    add_mapping_annotation(
        atoms_state.AtomsState.m_def,
        OUT_KEY,
        '.particle_states',
    )


class AtomsState(atoms_state.AtomsState):
    add_mapping_annotation(
        atoms_state.AtomsState.atomic_number,
        OUT_KEY,
        '.atomic_number',
    )

    add_mapping_annotation(
        atoms_state.AtomsState.chemical_symbol,
        OUT_KEY,
        '.chemical_symbol',
    )


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.value, OUT_KEY, '.value')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.value')


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.total_energies')
    add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.total_forces')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.functional_key, OUT_KEY, '.functional_key'
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.name, OUT_KEY, '.name')
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.xc')


class Simulation(general.Simulation):
    add_mapping_annotation(
        general.Simulation.program, OUT_KEY, ('get_program_data', ['.@'])
    )
    add_mapping_annotation(
        general.Simulation.model_system,
        OUT_KEY,
        ('get_atoms', ['.@']),
    )
    add_mapping_annotation(general.Simulation.outputs, OUT_KEY, ('get_outputs', ['.@']))
    add_mapping_annotation(
        model_method.DFT.m_def,
        OUT_KEY,
        ('get_dft_method', ['.@']),
    )


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')

try:
    m_package.__init_metainfo__()
except Exception:
    pass
