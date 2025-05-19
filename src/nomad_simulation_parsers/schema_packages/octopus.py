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
        dict(out=Mapper(mapper='.Version'))
    )


class XCFunctional(model_method.XCFunctional):
    model_method.XCFunctional.libxc_name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))


class DFT(model_method.DFT):
    model_method.DFT.xc_functionals.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_xc_functionals', ['.theory_level']))))


class ModelSystem(model_system.ModelSystem):
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.positions')))


class TotalEnergy(outputs.TotalEnergy):
    outputs.TotalEnergy.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.energy'), info=Mapper(mapper='.Total || .value')))


class TotalForce(outputs.TotalForce):
    outputs.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.forces')))


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    outputs.ElectronicEigenvalues.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            info=Mapper(mapper='.eigenvalues'),
            eigenvalues=Mapper(mapper='.eigenvalues'),
        )
    )
    outputs.ElectronicEigenvalues.occupation.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            info=Mapper(mapper='.ocupations'), eigenvalues=Mapper(mapper='.occupations')
        )
    )


class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@'), info=Mapper(mapper='.energies')))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.@')))
    outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            info=Mapper(mapper=('get_eigenvalues', ['eigenvalues.eigenvalues'])),
            eigenvalues=Mapper(mapper=('get_eigenvalues', ['eigenvalues'])),
        )
    )


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_header', ['.header']))))
    model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.@'))
    )
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_systems', ['.minimization']))))
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=('get_outputs', ['.minimization || .time_dependent[].iteration'])
            ),
            info=Mapper(mapper='.@'),
            eigenvalues=Mapper(mapper='.@'),
        )
    )


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        out=Mapper(mapper='@'),
        info=Mapper(mapper='@'),
        eigenvalues=Mapper(mapper='@'),
    )
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
