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
        dict(out=Mapper(mapper=('get_version', ['.program_name_version'])))
    )
    general.Program.datetime.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_datetime', ['.start_date_time']))))


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(out=Mapper(mapper='.XC_functional_name')))


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(out=Mapper(mapper=('get_xc_functionals', ['.xc_functional']))))


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.chemical_symbol.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))


class ModelSystem(model_system.ModelSystem):
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=('get_value', ['.@'], dict(key='labels_positions.positions'))
            )
        )
    )
    model_system.AtomsState.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=(
                    'get_value',
                    ['.@'],
                    dict(key='labels_positions.labels', units=''),
                )
            )
        )
    )
    model_system.AtomicCell.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(out=Mapper(mapper=('get_value', ['.@'], dict(key='simulation_cell'))))
    )


class TotalEmergy(outputs.TotalEnergy):
    outputs.TotalEnergy.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.value || .energy_total', unit='rydberg')))
    outputs.TotalEnergy.contributions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_energy_contributions', ['.@']))))


class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.energies')))


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.header')))
    model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.header'))
    )
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@', cache=True)))
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@', cache=True)))


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(out=Mapper(mapper='@'))
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
