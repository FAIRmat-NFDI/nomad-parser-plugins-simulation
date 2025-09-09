from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations

m_package = SchemaPackage()


class Program(general.Program):
    add_mapping_annotations(
        general.Program.version, 'out', ('get_version', ['.program_name_version'])
    )
    add_mapping_annotations(
        general.Program.datetime, 'out', ('get_datetime', ['.start_date_time'])
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotations(
        model_method.XCComponent.canonical_label, 'out', '.XC_functional_name'
    )


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotations(
        model_method.XCFunctional.components,
        'out',
        ('get_xc_functionals', ['.xc_functional']),
    )


class DFT(model_method.DFT):
    add_mapping_annotations(model_method.DFT.xc, 'out', '@')


class AtomsState(model_system.AtomsState):
    add_mapping_annotations(model_system.AtomsState.chemical_symbol, 'out', '.@')


class AtomicCell(model_system.AtomicCell):
    add_mapping_annotations(model_system.AtomicCell.lattice_vectors, 'out', '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotations(
        model_system.ModelSystem.positions,
        'out',
        ('get_value', ['.@'], dict(key='labels_positions.positions')),
    )
    add_mapping_annotations(
        model_system.AtomsState.m_def,
        'out',
        (
            'get_value',
            ['.@'],
            dict(key='labels_positions.labels', units=''),
        ),
    )
    add_mapping_annotations(
        model_system.AtomicCell.m_def,
        'out',
        ('get_value', ['.@'], dict(key='simulation_cell')),
    )


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotations(
        outputs.TotalEnergy.value, 'out', '.value || .energy_total', unit='rydberg'
    )
    add_mapping_annotations(
        outputs.TotalEnergy.contributions, 'out', ('get_energy_contributions', ['.@'])
    )


class Outputs(outputs.Outputs):
    add_mapping_annotations(outputs.Outputs.total_energies, 'out', '.energies')


class Simulation(general.Simulation):
    add_mapping_annotations(general.Simulation.program, 'out', '.header')
    add_mapping_annotations(model_method.DFT.m_def, 'out', '.header')
    add_mapping_annotations(general.Simulation.model_system, 'out', '.@')
    add_mapping_annotations(general.Simulation.outputs, 'out', '.@')


add_mapping_annotations(general.Simulation.m_def, 'out', '@')

try:
    m_package.__init_metainfo__()
except Exception:
    pass
