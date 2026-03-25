from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    force_field,
    general,
    model_system,
    outputs,
)
from nomad_simulations.schema_packages.workflow import (
    geometry_optimization,
    molecular_dynamics,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()


LOG_KEY = 'gromacs_log'
TPR_KEY = 'gromacs_tpr'
EDR_KEY = 'gromacs_edr'


# =============================================================================
# FORCE FIELD PARSING - TPR FILE EXTENSION ROADMAP
# =============================================================================
#
# CURRENT STATE:
# -------------
# - Bond topology (bond_list) is extracted from TPR via MDAnalysis
# - Force field parameters are read but not fully connected to topology
# - Contributions have particle_indices but lack parameter values
#
# TPR FILE CONTAINS:
# -----------------
# 1. Force field parameter sets (functypes array)
#    - Each set has type ID (F_BONDS, F_ANGLES, etc.) and parameter values
#    - Example: functype=F_BONDS, params=[k=500.0, r0=0.15]
#
# 2. Interaction lists (ilist section)
#    - Maps atom pairs/triplets/quadruplets to parameter set indices
#    - Example: bond 0-1 uses parameter set 3
#
# 3. Atom properties (charges, masses, atom types)
#
# FOR FULL IMPLEMENTATION:
# -----------------------
# 1. Extend MDAnalysis or implement custom TPR reader to access ilist
# 2. Create mapping: interaction -> parameter set -> actual values
# 3. Populate force_field.Potential subclasses:
#    - HarmonicBond: bond_constant, reference_bond_length
#    - MorseBond: dissociation_energy, width_parameter, reference_bond_length
#    - HarmonicAngle: angle_constant, reference_angle
#    - LennardJones: sigma, epsilon (or C6, C12)
#
# 4. Add schema annotations below to map TPR data to ForceField sections
#
# GROMACS TYPE -> NOMAD CLASS MAPPING:
# ------------------------------------
# F_BONDS (0) -> force_field.HarmonicBond
# F_G96BONDS (1) -> force_field.HarmonicBond (different units)
# F_MORSE (2) -> force_field.MorseBond
# F_ANGLES (10) -> force_field.HarmonicAngle
# F_PDIHS (19) -> force_field.ProperDihedral
# F_RBDIHS (27) -> force_field.RyckaertBellemansDihedral
# F_LJ (37) -> force_field.LennardJones
# F_LJ14 (45) -> force_field.LennardJones (1-4 interactions)
# F_SETTLE (64) -> force_field.Constraint (water geometry)
# F_CONSTR (62) -> force_field.Constraint
#
# Example annotations (for future implementation):
# ------------------------------------------------
# class ForceField(force_field.ForceField):
#     add_mapping_annotation(
#         force_field.ForceField.contributions,
#         TPR_KEY,
#         ('get_force_field_contributions_with_params', []),
#     )
#
# class HarmonicBond(force_field.HarmonicBond):
#     add_mapping_annotation(
#         force_field.HarmonicBond.bond_constant,
#         TPR_KEY,
#         '.bond_parameters.k'
#     )
#     add_mapping_annotation(
#         force_field.HarmonicBond.reference_bond_length,
#         TPR_KEY,
#         '.bond_parameters.r0'
#     )
# =============================================================================


class Program(general.Program):
    add_mapping_annotation(
        general.Program.version, LOG_KEY, ('get_version', ['.version'])
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.label, TPR_KEY, '.@')


add_mapping_annotation(
    model_system.ModelSystem.lattice_vectors, TPR_KEY, '.lattice_vectors'
)
add_mapping_annotation(
    model_system.ModelSystem.periodic_boundary_conditions, LOG_KEY, '.pbc'
)


class ModelSystem(model_system.ModelSystem):
    """
    GROMACS model system with topology information extracted from TPR files.

    IMPLEMENTED:
    - positions, velocities: Particle coordinates and velocities
    - bond_list: Bond topology (atom index pairs) extracted from MDAnalysis
    - particle_states (AtomsState): Atom labels and types
    - sub_systems: Hierarchical molecular structure (molecule groups → molecules →
      monomers)
    """

    add_mapping_annotation(model_system.ModelSystem.n_particles, TPR_KEY, '.n_atoms')
    add_mapping_annotation(model_system.ModelSystem.velocities, TPR_KEY, '.velocities')
    add_mapping_annotation(model_system.ModelSystem.positions, TPR_KEY, '.positions')
    add_mapping_annotation(model_system.ModelSystem.bond_list, TPR_KEY, '.bond_list')
    add_mapping_annotation(model_system.AtomsState.m_def, TPR_KEY, '.labels')


# ROOT annotation for nested ModelSystem instances
# When MappingParser creates ModelSystem from subsystems list,
# use that dict as source root
add_mapping_annotation(model_system.ModelSystem.m_def, TPR_KEY, '@')

# Subsystem hierarchy annotations (apply to all ModelSystem instances including
# subsystems)
# sub_systems: recursively extract from nested dicts via function call
# Pass '.@' as first argument (current node dict) to function
add_mapping_annotation(
    model_system.ModelSystem.sub_systems,
    TPR_KEY,
    ('get_subsystems_from_dict', ['.@']),
)
add_mapping_annotation(model_system.ModelSystem.name, TPR_KEY, '.name')
add_mapping_annotation(
    model_system.ModelSystem.composition_formula, TPR_KEY, '.composition_formula'
)
add_mapping_annotation(
    model_system.ModelSystem.particle_indices, TPR_KEY, '.particle_indices'
)
add_mapping_annotation(model_system.ModelSystem.branch_label, TPR_KEY, '.branch_label')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.name, LOG_KEY, '.label')
    add_mapping_annotation(outputs.TotalEnergy.name, EDR_KEY, '.label')
    add_mapping_annotation(outputs.TotalEnergy.value, LOG_KEY, '.value')
    add_mapping_annotation(outputs.TotalEnergy.value, EDR_KEY, '.value')
    add_mapping_annotation(outputs.TotalEnergy.contributions, LOG_KEY, '.contributions')
    add_mapping_annotation(outputs.TotalEnergy.contributions, EDR_KEY, '.contributions')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, TPR_KEY, '.@')


class Outpus(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, LOG_KEY, '.energy')
    add_mapping_annotation(outputs.Outputs.total_energies, EDR_KEY, '.energy')
    add_mapping_annotation(outputs.Outputs.model_system_ref, LOG_KEY, '.system_ref')
    add_mapping_annotation(outputs.Outputs.model_system_ref, EDR_KEY, '.system_ref')
    add_mapping_annotation(outputs.Outputs.total_forces, TPR_KEY, '.forces')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, LOG_KEY, '.header')
    add_mapping_annotation(
        general.Simulation.model_system, LOG_KEY, ('get_configurations', [])
    )
    add_mapping_annotation(
        general.Simulation.model_system, TPR_KEY, ('get_configurations', [])
    )
    add_mapping_annotation(general.Simulation.outputs, LOG_KEY, ('get_outputs', []))
    add_mapping_annotation(general.Simulation.outputs, TPR_KEY, ('get_outputs', []))
    add_mapping_annotation(general.Simulation.outputs, EDR_KEY, ('get_outputs', ['.@']))
    add_mapping_annotation(general.Simulation.model_method, TPR_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_method, LOG_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, LOG_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, TPR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EDR_KEY, '@')


class GeometryOptimizationModel(geometry_optimization.GeometryOptimization):
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationMethod.optimization_method,
        LOG_KEY,
        ('get_integrator_type', ['.input_parameters.integrator']),
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationMethod.n_steps_maximum,
        LOG_KEY,
        '.input_parameters.nsteps',
    )


class GeometryOptimizationResults(geometry_optimization.GeometryOptimizationResults):
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.energies,
        LOG_KEY,
        ('get_energies', []),
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.energies,
        EDR_KEY,
        ('get_energies', []),
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.final_force_maximum,
        LOG_KEY,
        'maximum_force',
        unit='kilojoule/avogadro_number/nanometer',
    )


class GeometryOptimization(geometry_optimization.GeometryOptimization):
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationMethod.m_def, LOG_KEY, '.@'
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.m_def, LOG_KEY, '.@'
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.m_def, EDR_KEY, '.@'
    )


class MolecularDynamicsModel(molecular_dynamics.MolecularDynamicsMethod):
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.integrator_type,
        LOG_KEY,
        ('get_integrator_type', ['.input_parameters.integrator']),
    )
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.integration_timestep,
        LOG_KEY,
        '.input_parameters.dt',
        unit='picosecond',
    )
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.n_steps,
        LOG_KEY,
        '.input_parameters.nsteps',
    )

    # Trajectory output frequencies
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.coordinate_save_frequency,
        LOG_KEY,
        ('get_coordinate_save_frequency', ['.input_parameters']),
    )
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.velocity_save_frequency,
        LOG_KEY,
        '.input_parameters.nstvout',
    )
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.force_save_frequency,
        LOG_KEY,
        '.input_parameters.nstfout',
    )
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.thermodynamics_save_frequency,
        LOG_KEY,
        '.input_parameters.nstenergy',
    )

    # Thermodynamic ensemble
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.thermodynamic_ensemble,
        LOG_KEY,
        ('get_thermodynamic_ensemble', ['.input_parameters']),
    )

    # Thermostat subsection
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.thermostat_parameters,
        LOG_KEY,
        '@',
    )

    # Barostat subsection
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.barostat_parameters,
        LOG_KEY,
        '@',
    )


## ThermostatParameters annotations

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.thermostat_type,
    LOG_KEY,
    ('get_thermostat_type', ['.input_parameters.tcoupl']),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.reference_temperature,
    LOG_KEY,
    ('get_reference_temperature', ['.input_parameters']),
    unit='kelvin',
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.coupling_constant,
    LOG_KEY,
    ('get_thermostat_coupling_constant', ['.input_parameters']),
    unit='picosecond',
)

## BarostatParameters annotations

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.barostat_type,
    LOG_KEY,
    ('get_barostat_type', ['.input_parameters.pcoupl']),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.coupling_type,
    LOG_KEY,
    ('get_barostat_coupling_type', ['.input_parameters.pcoupltype']),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.reference_pressure,
    LOG_KEY,
    ('get_reference_pressure', ['.input_parameters']),
    unit='bar',
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.coupling_constant,
    LOG_KEY,
    ('get_barostat_coupling_constant', ['.input_parameters']),
    unit='picosecond',
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.compressibility,
    LOG_KEY,
    ('get_compressibility', ['.input_parameters']),
    unit='1/bar',
)


class MolecularDynamicsResults(molecular_dynamics.MolecularDynamicsResults):
    # parse from xvg
    pass


class MolecularDynamics(molecular_dynamics.MolecularDynamics):
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsMethod.m_def, LOG_KEY, '.@'
    )


# Workflow
add_mapping_annotation(geometry_optimization.GeometryOptimization.m_def, LOG_KEY, '@')
add_mapping_annotation(geometry_optimization.GeometryOptimization.m_def, EDR_KEY, '@')
add_mapping_annotation(molecular_dynamics.MolecularDynamics.m_def, LOG_KEY, '@')


# Force Field
class ForceField(force_field.ForceField):
    add_mapping_annotation(
        force_field.ForceField.contributions,
        TPR_KEY,
        ('get_force_field_contributions', []),
    )


class ForceCalculations(force_field.ForceCalculations):
    add_mapping_annotation(
        force_field.ForceCalculations.vdw_cutoff,
        LOG_KEY,
        '.input_parameters.rvdw',
        unit='nanometer',
    )
    add_mapping_annotation(
        force_field.ForceCalculations.coulomb_cutoff,
        LOG_KEY,
        '.input_parameters.rcoulomb',
        unit='nanometer',
    )
    add_mapping_annotation(
        force_field.ForceCalculations.coulomb_type,
        LOG_KEY,
        ('get_coulomb_type', ['.input_parameters.coulombtype']),
    )
    add_mapping_annotation(
        force_field.ForceCalculations.neighbor_update_frequency,
        LOG_KEY,
        '.input_parameters.nstlist',
    )


# add_mapping_annotation(force_field.ForceField.m_def, TPR_KEY, '@')
add_mapping_annotation(force_field.Potential.m_def, TPR_KEY, '@')
add_mapping_annotation(force_field.ForceCalculations.m_def, LOG_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
