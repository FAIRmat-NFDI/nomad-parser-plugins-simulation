from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()


OUT_KEY = 'crystal_out'
F25_KEY = 'crystal_f25'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.label')
    add_mapping_annotation(model_system.AtomsState.atomic_number, OUT_KEY, '.number')


class Representation(model_system.Representation):
    # TODO the or || operator does not seem to work
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, OUT_KEY, '.lattice_vectors'
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, OUT_KEY, '.@')
    add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, OUT_KEY, '.atoms')


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, OUT_KEY, '.name')


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         out=Mapper(mapper='.name')
#     )


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         out=Mapper(mapper=('get_xc_functionals', ['.@']))
#     )


# class TotalEnergy(outputs.TotalEnergy):
#     outputs.TotalEnergy.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         out=Mapper(mapper='.energy')
#     )


class TotalForces(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.forces')


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(outputs.ElectronicDensityOfStates.value, F25_KEY, '.values')


class Outputs(outputs.Outputs):
    # outputs.Outputs.total_energies.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    #     out=Mapper(mapper='.@')
    # )
    outputs.Outputs.total_forces.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.@')
    )
    outputs.Outputs.electronic_dos.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        f25=Mapper(mapper=('get_dos', ['.dos']))
    )
    outputs.Outputs.electronic_band_structures.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(out=Mapper(mapper=('get_band_structures', ['.band_structure'])))
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.wall_start, OUT_KEY, ('to_unix_time', ['.start_timestamp'])
    )
    add_mapping_annotation(
        general.Simulation.wall_end, OUT_KEY, ('to_unix_time', ['.end_timestamp'])
    )
    add_mapping_annotation(
        general.Simulation.model_system, OUT_KEY, ('get_systems', ['.@'])
    )
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.dft')
    add_mapping_annotation(general.Simulation.outputs, OUT_KEY, ('get_outputs', ['.@']))
    add_mapping_annotation(general.Simulation.outputs, F25_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, F25_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
