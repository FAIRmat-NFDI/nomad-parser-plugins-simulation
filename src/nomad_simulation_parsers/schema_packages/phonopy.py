import numpy as np
from nomad.metainfo import Quantity, SchemaPackage, SubSection
from nomad_simulations.schema_packages import general, model_system, variables
from nomad_simulations.schema_packages.model_method import ModelMethod
from nomad_simulations.schema_packages.outputs import Outputs
from nomad_simulations.schema_packages.physical_property import PhysicalProperty
from nomad_simulations.schema_packages.properties import (
    DOSProfile,
    HeatCapacity,
    HelmholtzFreeEnergy,
)
from nomad_simulations.schema_packages.variables import (
    Temperature,
    Variables,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()
PHONOPY_KEY = 'phonopy_yaml'


class PhonopyMethod(ModelMethod):
    """Parameters of the finite-displacement lattice-dynamics model."""

    displacement = Quantity(
        type=np.float64,
        unit='meter',
        description='Amplitude of the atomic displacement used for the supercell.',
    )
    add_mapping_annotation(displacement, PHONOPY_KEY, '.displacement')

    symmetry_tolerance = Quantity(
        type=np.float64,
        unit='meter',
        description='Cartesian tolerance used to determine crystal symmetry.',
    )
    add_mapping_annotation(symmetry_tolerance, PHONOPY_KEY, '.symmetry_tolerance')

    force_calculator = Quantity(
        type=str,
        description='Program used to calculate the displaced-supercell forces.',
    )
    add_mapping_annotation(force_calculator, PHONOPY_KEY, '.force_calculator')

    mesh_density = Quantity(
        type=np.float64,
        unit='1 / meter ** 3',
        description='Density of the reciprocal-space mesh used by Phonopy.',
    )
    add_mapping_annotation(mesh_density, PHONOPY_KEY, '.mesh_density')

    with_non_analytic_correction = Quantity(
        type=bool,
        description='Whether non-analytical term corrections are applied.',
    )
    add_mapping_annotation(
        with_non_analytic_correction,
        PHONOPY_KEY,
        '.with_non_analytic_correction',
    )


add_mapping_annotation(PhonopyMethod.name, PHONOPY_KEY, '.name')
add_mapping_annotation(PhonopyMethod.type, PHONOPY_KEY, '.type')
add_mapping_annotation(PhonopyMethod.contributions, PHONOPY_KEY, '.contributions')


class PhonopyProgram(general.Program):
    add_mapping_annotation(general.Program.name, PHONOPY_KEY, '.name')
    add_mapping_annotation(general.Program.version, PHONOPY_KEY, '.version')


class PhonopyAtomsState(model_system.AtomsState):
    add_mapping_annotation(
        model_system.AtomsState.chemical_symbol, PHONOPY_KEY, '.chemical_symbol'
    )


class PhonopyRepresentation(model_system.AlternativeRepresentation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, PHONOPY_KEY, '.lattice_vectors'
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        PHONOPY_KEY,
        '.periodic_boundary_conditions',
    )
    add_mapping_annotation(
        model_system.AlternativeRepresentation.supercell_matrix,
        PHONOPY_KEY,
        '.supercell_matrix',
    )


class PhonopyModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions, PHONOPY_KEY, '.positions'
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def, PHONOPY_KEY, '.particle_states'
    )
    add_mapping_annotation(
        model_system.AlternativeRepresentation.m_def, PHONOPY_KEY, '.representations'
    )


class ForceConstants(PhysicalProperty):
    """Second derivatives with respect to pairs of atomic Cartesian coordinates."""

    value = Quantity(
        type=np.float64,
        unit='joule / meter ** 2',
        shape=['*', '*', 3, 3],
    )
    add_mapping_annotation(value, PHONOPY_KEY, '.value')


class QPath(Variables):
    """Reciprocal-space path points in fractional coordinates."""

    points = Quantity(type=np.float64, shape=['n_points', 3])
    add_mapping_annotation(points, PHONOPY_KEY, '.points')


class PhononBandStructure(PhysicalProperty):
    """Phonon mode energies along one reciprocal-space path segment."""

    n_bands = Quantity(type=np.int32)
    value = Quantity(
        type=np.float64,
        unit='joule',
        shape=['*', 'n_bands'],
    )
    endpoints_labels = Quantity(type=str, shape=[2])
    q_path = SubSection(sub_section=QPath.m_def)
    add_mapping_annotation(value, PHONOPY_KEY, '.value')
    add_mapping_annotation(n_bands, PHONOPY_KEY, '.n_bands')
    add_mapping_annotation(endpoints_labels, PHONOPY_KEY, '.endpoints_labels')
    add_mapping_annotation(q_path, PHONOPY_KEY, '.q_path')


class PhononDensityOfStates(DOSProfile):
    """Phonon density of states as a function of phonon mode energy."""

    value = Quantity(type=np.float64, unit='1 / joule', shape=['*'])
    add_mapping_annotation(value, PHONOPY_KEY, '.value')


add_mapping_annotation(variables.Energy2.points, PHONOPY_KEY, '.frequencies.points')


class VibrationalFreeEnergy(HelmholtzFreeEnergy):
    """Constant-volume vibrational free energy over a temperature grid."""

    value = Quantity(type=np.float64, unit='joule', shape=['*'])
    add_mapping_annotation(value, PHONOPY_KEY, '.value')
    temperatures = SubSection(sub_section=Temperature.m_def)


add_mapping_annotation(
    variables.Temperature.points, PHONOPY_KEY, '.temperatures.points'
)


class VibrationalHeatCapacity(HeatCapacity):
    """Constant-volume vibrational heat capacity over a temperature grid."""

    value = Quantity(type=np.float64, unit='joule / kelvin', shape=['*'])
    add_mapping_annotation(value, PHONOPY_KEY, '.value')
    temperatures = SubSection(sub_section=Temperature.m_def)


add_mapping_annotation(
    variables.Temperature.points, PHONOPY_KEY, '.temperatures.points'
)


class PhonopyOutputs(Outputs):
    force_constants = SubSection(sub_section=ForceConstants.m_def, repeats=True)
    add_mapping_annotation(force_constants, PHONOPY_KEY, ('get_force_constants', []))
    phonon_band_structures = SubSection(
        sub_section=PhononBandStructure.m_def, repeats=True
    )
    add_mapping_annotation(
        phonon_band_structures, PHONOPY_KEY, ('get_phonon_band_structures', [])
    )
    phonon_dos = SubSection(sub_section=PhononDensityOfStates.m_def, repeats=True)
    add_mapping_annotation(phonon_dos, PHONOPY_KEY, ('get_phonon_dos', []))
    vibrational_free_energies = SubSection(
        sub_section=VibrationalFreeEnergy.m_def, repeats=True
    )
    add_mapping_annotation(
        vibrational_free_energies,
        PHONOPY_KEY,
        ('get_vibrational_thermodynamics', [], {'key': 'vibrational_free_energies'}),
        cache=True,
    )
    vibrational_heat_capacities = SubSection(
        sub_section=VibrationalHeatCapacity.m_def, repeats=True
    )
    add_mapping_annotation(
        vibrational_heat_capacities,
        PHONOPY_KEY,
        (
            'get_vibrational_thermodynamics',
            [],
            {'key': 'vibrational_heat_capacities'},
        ),
        cache=True,
    )

    n_imaginary_frequencies = Quantity(
        type=np.int32,
        description='Number of negative frequencies on the DOS sampling mesh.',
    )
    add_mapping_annotation(
        n_imaginary_frequencies, PHONOPY_KEY, ('get_imaginary_frequencies', [])
    )


class PhonopySimulation(general.Simulation):
    add_mapping_annotation(PhonopyMethod.m_def, PHONOPY_KEY, ('get_model_method', []))
    add_mapping_annotation(PhonopyOutputs.m_def, PHONOPY_KEY, '.@')


add_mapping_annotation(PhonopySimulation.program, PHONOPY_KEY, ('get_program', []))
add_mapping_annotation(
    PhonopySimulation.model_system, PHONOPY_KEY, ('get_model_system', [])
)
add_mapping_annotation(PhonopySimulation.m_def, PHONOPY_KEY, '@')


m_package.__init_metainfo__()
