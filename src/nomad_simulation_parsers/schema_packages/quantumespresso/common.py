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
XML_KEY = 'xml'


class Program(general.Program):
    add_mapping_annotations(
        general.Program.version, OUT_KEY, ('get_version', ['.program_name_version'])
    )
    add_mapping_annotations(general.Program.version, XML_KEY, '.creator."@VERSION"')
    add_mapping_annotations(
        general.Program.datetime, OUT_KEY, ('get_datetime', ['.start_date_time'])
    )
    add_mapping_annotations(
        general.Program.datetime,
        XML_KEY,
        ('get_datetime', ['.created."@DATE"', '.created."@TIME"']),
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
    add_mapping_annotations(model_system.AtomsState.chemical_symbol, OUT_KEY, '.@')
    add_mapping_annotations(
        model_system.AtomsState.chemical_symbol, XML_KEY, '."@name"'
    )


class Representation(model_system.Representation):
    add_mapping_annotations(model_system.Representation.lattice_vectors, 'out', '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotations(
        model_system.ModelSystem.positions,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='labels_positions.positions')),
    )
    add_mapping_annotations(
        model_system.ModelSystem.positions,
        XML_KEY,
        (
            'apply_unit',
            ['.atomic_structure.atomic_positions.atom[].__value'],
            dict(name='length'),
        ),
    )
    add_mapping_annotations(
        model_system.AtomsState.m_def,
        OUT_KEY,
        (
            'get_value',
            ['.@'],
            dict(key='labels_positions.labels', units=''),
        ),
    )
    add_mapping_annotations(
        model_system.Representation.m_def,
        'out',
        ('get_value', ['.@'], dict(key='simulation_cell')),
    )
    add_mapping_annotations(
        model_system.AtomicCell.m_def, XML_KEY, '.atomic_structure.cell'
    )


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotations(
        outputs.TotalEnergy.value, OUT_KEY, '.value || .energy_total', unit='rydberg'
    )
    add_mapping_annotations(outputs.TotalEnergy.value, XML_KEY, '.value || .etot')
    add_mapping_annotations(
        outputs.TotalEnergy.contributions, OUT_KEY, ('get_energy_contributions', ['.@'])
    )
    add_mapping_annotations(
        outputs.TotalEnergy.contributions, XML_KEY, ('get_energy_contributions', ['.@'])
    )
    add_mapping_annotations(outputs.TotalEnergy.name, XML_KEY, '.name')


class Outputs(outputs.Outputs):
    add_mapping_annotations(outputs.Outputs.total_energies, OUT_KEY, '.energies')
    add_mapping_annotations(outputs.Outputs.total_energies, XML_KEY, '.total_energy')


class Simulation(general.Simulation):
    add_mapping_annotations(general.Simulation.program, OUT_KEY, '.header')
    add_mapping_annotations(general.Simulation.program, XML_KEY, '.general_info')
    add_mapping_annotations(model_method.DFT.m_def, OUT_KEY, '.header')
    add_mapping_annotations(model_method.DFT.m_def, XML_KEY, '.input')
    add_mapping_annotations(general.Simulation.model_system, OUT_KEY, '.@')
    add_mapping_annotations(general.Simulation.outputs, OUT_KEY, '.@')


add_mapping_annotations(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotations(general.Simulation.m_def, XML_KEY, '@')

try:
    m_package.__init_metainfo__()
except Exception:
    pass
