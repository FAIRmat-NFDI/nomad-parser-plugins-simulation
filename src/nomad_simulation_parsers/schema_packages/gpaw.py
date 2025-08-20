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
    general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(gpw=Mapper(mapper='.program_version'))
    )


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.unitcell')))
    model_system.AtomicCell.periodic_boundary_conditions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.boundary_conditions')))


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.chemical_symbol.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.@')))


class ModelSystem(model_system.ModelSystem):
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.atom_positions')))
    model_system.AtomicCell.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.@')))
    model_system.AtomsState.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.labels')))


class XCFunctional(model_method.XCFunctional):
    model_method.XCFunctional.libxc_name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.@')))


class DFT(model_method.DFT):
    model_method.DFT.xc_functionals.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.xcfunctional')))


class TotalEnergy(outputs.TotalEnergy):
    outputs.TotalEnergy.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.total || .value')))
    outputs.TotalEnergy.name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.name')))
    outputs.TotalEnergy.contributions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.contributions')))


class TotalForce(outputs.TotalForce):
    outputs.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.total || .value')))


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    outputs.ElectronicEigenvalues.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.eigenvalues')))
    outputs.ElectronicEigenvalues.occupation.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.occupations')))


class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper=('get_energies', []))))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper=('get_forces', []))))
    outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper=('get_eigenvalues', []))))


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.@')))
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.@')))
    model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(gpw=Mapper(mapper='.@'))
    )
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(gpw=Mapper(mapper='.@')))


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(gpw=Mapper(mapper='@'))
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
