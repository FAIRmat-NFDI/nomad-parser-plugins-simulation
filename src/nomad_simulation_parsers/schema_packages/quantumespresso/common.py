from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'quantumespresso_out'


class Program(general.Program):
    add_mapping_annotation(
        general.Program.version, OUT_KEY, ('get_version', ['.program_name_version'])
    )
    add_mapping_annotation(
        general.Program.datetime, OUT_KEY, ('get_datetime', ['.start_date_time'])
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(
        model_method.XCComponent.canonical_label, OUT_KEY, '.XC_functional_name'
    )


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        OUT_KEY,
        ('get_xc_functionals', ['.xc_functional']),
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.@')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.@')


class Representation(model_system.Representation):
    add_mapping_annotation(model_system.Representation.lattice_vectors, OUT_KEY, '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='labels_positions.positions')),
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='labels_positions.labels')),
    )
    add_mapping_annotation(
        model_system.Representation.m_def,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='simulation_cell')),
    )


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(
        outputs.TotalEnergy.value, OUT_KEY, '.value || .energy_total', unit='rydberg'
    )
    add_mapping_annotation(
        outputs.TotalEnergy.contributions, OUT_KEY, ('get_energy_contributions', ['.@'])
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.energies')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.header')
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.header')
    add_mapping_annotation(general.Simulation.model_system, OUT_KEY, '.@', cache=True)
    add_mapping_annotation(general.Simulation.outputs, OUT_KEY, '.@', cache=True)


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '.@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
