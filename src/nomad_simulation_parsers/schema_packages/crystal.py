from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
    variables,
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
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        OUT_KEY,
        ('get_periodic_boundary_conditions', ['.lattice_vectors']),
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, OUT_KEY, '.@')
    add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, OUT_KEY, '.atoms')


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, OUT_KEY, '.name')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components, OUT_KEY, ('get_xc_functionals', ['.dft'])
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.@')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.value, OUT_KEY, '.energy')


class TotalForces(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, OUT_KEY, '.forces')


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(outputs.ElectronicDensityOfStates.value, F25_KEY, '.values')


class Energy2(variables.Energy2):
    add_mapping_annotation(variables.Energy2.points, F25_KEY, '.energies')


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    add_mapping_annotation(
        outputs.ElectronicBandStructure.value,
        F25_KEY,
        '.value',
        unit='hartree',
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.scf_steps, OUT_KEY, '.scf_steps')
    # TODO(legacy-parity): Electronic DOS and band structures are currently not
    # functional for the Crystal fixture set that only provides `si.f9`/`si.f98`.
    # Current mapping expects fort.25-compatible payload (`.dos`/`.segments`).
    # Re-enable these mappings once parser-side fallback for f9/f98 is implemented.
    # add_mapping_annotation(
    #     outputs.Outputs.electronic_dos, F25_KEY, ('get_dos', ['.dos'])
    # )
    # add_mapping_annotation(
    #     outputs.Outputs.electronic_band_structures,
    #     F25_KEY,
    #     ('get_band_structures', ['.segments']),
    # )
    # TODO(legacy-parity): no explicit crystal band-gap section was populated in
    # legacy parser output; keep unmapped until dedicated legacy-equivalent source
    # is confirmed.


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


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(outputs.SCFSteps.energies_total, OUT_KEY, '.energies_total')
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, OUT_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, OUT_KEY, '.code_specific_quantities'
    )


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, F25_KEY, '@')
add_mapping_annotation(variables.Energy2.m_def, F25_KEY, '.@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
