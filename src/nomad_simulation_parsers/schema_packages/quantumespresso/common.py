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
XML_KEY = 'quantumespresso_xml'
GIPAW_OUT_KEY = 'quantumespresso_gipaw_out'
GIPAW_XML_KEY = 'quantumespresso_gipaw_xml'
DOS_KEY = 'quantumespresso_dos'


class Program(general.Program):
    add_mapping_annotation(
        general.Program.version, OUT_KEY, ('get_version', ['.program_name_version'])
    )
    add_mapping_annotation(general.Program.version, XML_KEY, '.creator."@VERSION"')
    add_mapping_annotation(
        general.Program.datetime, OUT_KEY, ('get_datetime', ['.start_date_time'])
    )
    add_mapping_annotation(
        general.Program.datetime,
        XML_KEY,
        ('get_datetime', ['.created."@DATE"', '.created."@TIME"']),
    )


class DFT(model_method.DFT):
    # Materialize the `xc` subsection so its child `functional_key` mapper runs.
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.@')
    add_mapping_annotation(model_method.DFT.xc, XML_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    # The output prints a functional name or the four DFT slot codes; the XML
    # carries a clean functional name. Both resolve to a standard name that the
    # schema expands into components (family/kind) and `jacobs_ladder`.
    add_mapping_annotation(
        model_method.XCFunctional.functional_key,
        OUT_KEY,
        ('get_functional_key', ['.xc_functional']),
    )
    add_mapping_annotation(
        model_method.XCFunctional.functional_key, XML_KEY, '.dft.functional'
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.@')
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, XML_KEY, '."@name"')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='simulation_cell')),
    )
    add_mapping_annotation(
        model_system.Representation.lattice_vectors,
        XML_KEY,
        ('apply_unit', ['.__value'], dict(name='length')),
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        OUT_KEY,
        ('get_periodic_boundary_conditions', ['.@']),
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        XML_KEY,
        ('get_periodic_boundary_conditions', ['.atomic_structure.cell']),
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='labels_positions.positions')),
    )
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        XML_KEY,
        (
            'apply_unit',
            ['.atomic_structure.atomic_positions.atom[].__value'],
            dict(name='length'),
        ),
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def,
        OUT_KEY,
        (
            'get_value',
            ['.@'],
            dict(key='labels_positions.labels', units=''),
        ),
    )
    add_mapping_annotation(
        model_system.Representation.m_def,
        OUT_KEY,
        ('get_value', ['.@'], dict(key='simulation_cell')),
    )
    add_mapping_annotation(
        model_system.Representation.m_def, XML_KEY, '.atomic_structure.cell'
    )


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(
        outputs.TotalEnergy.value, OUT_KEY, '.value || .energy_total', unit='rydberg'
    )
    add_mapping_annotation(outputs.TotalEnergy.value, XML_KEY, '.value || .etot')
    add_mapping_annotation(
        outputs.TotalEnergy.contributions, OUT_KEY, ('get_energy_contributions', ['.@'])
    )
    add_mapping_annotation(
        outputs.TotalEnergy.contributions, XML_KEY, ('get_energy_contributions', ['.@'])
    )
    add_mapping_annotation(outputs.TotalEnergy.name, XML_KEY, '.name')


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.energies')
    add_mapping_annotation(outputs.Outputs.total_energies, XML_KEY, '.total_energy')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.header')
    add_mapping_annotation(general.Simulation.program, XML_KEY, '.general_info')
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.header')
    add_mapping_annotation(model_method.DFT.m_def, XML_KEY, '.input')
    add_mapping_annotation(general.Simulation.model_system, OUT_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, OUT_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOS_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
