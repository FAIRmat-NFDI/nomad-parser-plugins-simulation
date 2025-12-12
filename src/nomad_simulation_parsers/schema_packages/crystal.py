from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

m_package = SchemaPackage()


class Program(general.Program):
    general.Program.version.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.program_version')
    )


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.chemical_symbol.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(out=Mapper(mapper='.label'))
    )
    model_system.AtomsState.atomic_number.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.number')
    )


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(
            # TODO the or || operator does not seem to work
            out=Mapper(mapper='.lattice_vectors')
        )
    )


class ModelSystem(model_system.ModelSystem):
    model_system.AtomicCell.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.@')
    )
    model_system.ModelSystem.positions.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.positions')
    )
    model_system.AtomsState.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.atoms')
    )


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
    outputs.TotalForce.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.forces')
    )


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    outputs.ElectronicDensityOfStates.value.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(f25=Mapper(mapper='.values'))
    )


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
    general.Simulation.program.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.@')
    )
    general.Simulation.wall_start.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper=('to_unix_time', ['.start_timestamp']))
    )
    general.Simulation.wall_end.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper=('to_unix_time', ['.end_timestamp']))
    )
    general.Simulation.model_system.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper=('get_systems', ['.@']))
    )
    model_method.DFT.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper='.dft')
    )
    general.Simulation.outputs.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        out=Mapper(mapper=('get_outputs', ['.@'])), f25=Mapper(mapper='.@')
    )


general.Simulation.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    out=Mapper(mapper='@'), f25=Mapper(mapper='@')
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
