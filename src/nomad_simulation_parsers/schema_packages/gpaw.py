from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

GPW_KEY = 'gpaw_gpw'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, GPW_KEY, '.program_version')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, GPW_KEY, '.unitcell'
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        GPW_KEY,
        '.boundary_conditions',
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, GPW_KEY, '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions, GPW_KEY, '.atom_positions'
    )
    add_mapping_annotation(model_system.Representation.m_def, GPW_KEY, '.@')
    add_mapping_annotation(model_system.AtomsState.m_def, GPW_KEY, '.labels')


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, GPW_KEY, '.@')


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(gpw=Mapper(mapper='.@')))


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(gpw=Mapper(mapper='.xcfunctional')))


# class TotalEnergy(outputs.TotalEnergy):
#     outputs.TotalEnergy.value.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(gpw=Mapper(mapper='.total || .value')))
#     outputs.TotalEnergy.name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(gpw=Mapper(mapper='.name')))
#     outputs.TotalEnergy.contributions.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(gpw=Mapper(mapper='.contributions')))


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, GPW_KEY, '.total || .value')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, GPW_KEY, '.eigenvalues')
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, GPW_KEY, '.occupations'
    )


class Outputs(outputs.Outputs):
    # outputs.Outputs.total_energies.m_annotations.setdefault(
    #     MAPPING_ANNOTATION_KEY, {}
    # ).update(dict(gpw=Mapper(mapper=('get_energies', []))))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper=('get_forces', []))))
    outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper=('get_eigenvalues', []))))


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
