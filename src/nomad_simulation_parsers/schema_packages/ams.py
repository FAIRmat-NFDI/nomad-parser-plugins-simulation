from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad_file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
    variables,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'ams_out'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.@')


class Representation(model_system.Representation):
    add_mapping_annotation(model_system.Representation.lattice_vectors, OUT_KEY, '.@')
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        OUT_KEY,
        ('get_periodic_boundary_conditions', ['.@']),
    )


class ModelSystem(model_system.ModelSystem):
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper='.labels_positions[1]')))
    model_system.AtomsState.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper='.labels_positions[0]')))
    model_system.Representation.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper='.lattice_vectors')))


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(out=Mapper(mapper='.@')))


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(
#         dict(
#             out=Mapper(
#                 mapper=('get_xc_functionals', ['.model_parameters.dft_potential'])
#             )
#         )
#     )


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(
        outputs.TotalEnergy.value, OUT_KEY, '.value || .energy_total'
    )
    add_mapping_annotation(
        outputs.TotalEnergy.contributions,
        OUT_KEY,
        ('get_contributions', ['.energies']),
    )


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.value || .forces_total')
    add_mapping_annotation(
        outputs.TotalForce.contributions, OUT_KEY, ('get_contributions', ['.forces'])
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, OUT_KEY, '.eigenvalues')
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, OUT_KEY, '.occupations'
    )

    # class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper='.@')))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper='.@')))
    outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            ams_out=Mapper(
                mapper=('get_eigenvalues', ['.eigenvalues || .band_energy_ranges'])
            )
        )
    )
    outputs.Outputs.scf_steps.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper=('get_scf_steps', ['.@']))))
    outputs.Outputs.electronic_band_gaps.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            ams_out=Mapper(
                mapper=(
                    'get_band_gaps',
                    ['.band_gap || .band_gap_info || .band_energy_ranges'],
                )
            )
        )
    )
    outputs.Outputs.electronic_dos.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(ams_out=Mapper(mapper=('get_dos', ['.total_dos']))))


class ElectronicBandGap(outputs.ElectronicBandGap):
    add_mapping_annotation(outputs.ElectronicBandGap.value, OUT_KEY, '.value')
    add_mapping_annotation(
        outputs.ElectronicBandGap.spin_channel, OUT_KEY, '.spin_channel'
    )


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value, OUT_KEY, '.value', unit='1 / hartree'
    )
    add_mapping_annotation(variables.Energy2.m_def, OUT_KEY, '.@')


class Energy2(variables.Energy2):
    add_mapping_annotation(
        variables.Energy2.points, OUT_KEY, '.energies', unit='hartree'
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        model_method.DFT.m_def,
        OUT_KEY,
        '.geometry_optimization || molecular_dynamics || .single_point',
    )
    add_mapping_annotation(
        general.Simulation.model_system,
        OUT_KEY,
        '.geometry_optimization.step|| molecular_dynamics.step || .single_point',
    )
    add_mapping_annotation(
        general.Simulation.outputs,
        OUT_KEY,
        '.geometry_optimization.step|| molecular_dynamics.step || .single_point',
    )


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, OUT_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, OUT_KEY, '.code_specific_quantities'
    )


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
