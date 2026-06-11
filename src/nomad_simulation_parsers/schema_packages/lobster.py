import numpy as np
from nomad.metainfo import (
    MSection,
    Quantity,
    SchemaPackage,
    Section,
    SectionProxy,
    SubSection,
)
from nomad_simulations.schema_packages import (
    basis_set,
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'lobster_out'
STRUCTURE_KEY = 'lobster_structure'
ICOXPLIST_KEY = 'lobster_icoxplist'


class x_lobster_section_cohp(MSection):
    """
    This is a section containing the crystal orbital hamilton population (COHP)
    and integrated COHP (iCOHP) values.
    """

    m_def = Section(validate=False)

    x_lobster_number_of_cohp_pairs = Quantity(
        type=int,
        description="""
        Number of atom pairs for which are the COHPs and iCOHPs calculated.
        """,
    )

    x_lobster_cohp_atom1_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Species and indices of the first atom for which is the specific COHP/iCOHP calculated
        """,
    )

    x_lobster_cohp_atom2_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Species and indices of the second atom for which is the specific COHP/iCOHP calculated
        """,
    )

    x_lobster_cohp_distances = Quantity(
        type=np.float64,
        unit='meter',
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Distance between atoms of the pair for which is the specific COHP/iCOHP calculated.
        """,
    )

    x_lobster_cohp_translations = Quantity(
        type=np.int32,
        shape=['x_lobster_number_of_cohp_pairs', 3],
        description="""
        Vector connecting the unit-cell of the first atom with the one of the second atom

        This is only used with LOBSTER versions 3.0.0 and above, older versions use
        x_lobster_cohp_number_of_bonds instead.
        """,
    )

    x_lobster_integrated_cohp_at_fermi_level = Quantity(
        type=np.float64,
        unit='joule',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_pairs'],
        description="""
        Calculated iCOHP values integrated up to the Fermi level.
        """,
    )

    x_lobster_number_of_cohp_values = Quantity(
        type=int,
        description="""
        Number of energy values for the COHP and iCOHP.
        """,
    )

    x_lobster_cohp_energies = Quantity(
        type=np.float64,
        unit='joule',
        shape=['x_lobster_number_of_cohp_values'],
        description="""
        Array containing the set of discrete energy values for COHP and iCOHP.
        """,
    )

    x_lobster_cohp_values = Quantity(
        type=np.float64,
        shape=[
            'x_lobster_number_of_cohp_pairs',
            'number_of_spin_channels',
            'x_lobster_number_of_cohp_values',
        ],
        description="""
        Calculated COHP values.
        """,
    )

    x_lobster_integrated_cohp_values = Quantity(
        type=np.float64,
        unit='joule',
        shape=[
            'x_lobster_number_of_cohp_pairs',
            'number_of_spin_channels',
            'x_lobster_number_of_cohp_values',
        ],
        description="""
        Calculated iCOHP values.
        """,
    )

    x_lobster_average_cohp_values = Quantity(
        type=np.float64,
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_values'],
        description="""
        Calculated COHP values averaged over all pairs.
        """,
    )

    x_lobster_average_integrated_cohp_values = Quantity(
        type=np.float64,
        unit='joule',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_values'],
        description="""
        Calculated iCOHP values averaged over all pairs.
        """,
    )

    x_lobster_cohp_number_of_bonds = Quantity(
        type=int,
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Number of bonds between first atom and the second atom (including
        the periodic images).

        This is only used in older LOBSTER versions, new versions print one line
        for every neighbor, so a pair which had x_lobster_icohp_number_of_bonds = 4
        in the old version would actually show as 4 lines in the ICOHPLIST or 4 columns
        in the COPHCAR in the new format.
        """,
    )
    x_lobster_cohp_orbital_per_label = SubSection(
        sub_section=SectionProxy('x_lobster_section_cohp_orbital_label'), repeats=True
    )


class Outputs(outputs.Outputs):
    x_lobster_section_cohp = SubSection(
        sub_section=SectionProxy('x_lobster_section_cohp')
    )
    add_mapping_annotation(x_lobster_section_cohp, ICOXPLIST_KEY, '.@')


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(
        model_system.AtomsState.atomic_number, STRUCTURE_KEY, '.number'
    )
    add_mapping_annotation(
        model_system.AtomsState.chemical_symbol, STRUCTURE_KEY, '.symbol'
    )


class ModelSystem(general.ModelSystem):
    add_mapping_annotation(
        general.ModelSystem.positions, STRUCTURE_KEY, '.positions', unit='angstrom'
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def, STRUCTURE_KEY, ('get_atoms', [])
    )
    add_mapping_annotation(
        general.ModelSystem.lattice_vectors, STRUCTURE_KEY, '.cell', unit='angstrom'
    )
    add_mapping_annotation(
        general.ModelSystem.periodic_boundary_conditions, STRUCTURE_KEY, '.pbc'
    )


class AtomCenteredBasisSet(basis_set.AtomCenteredBasisSet):
    add_mapping_annotation(
        basis_set.AtomCenteredBasisSet.basis_set, OUT_KEY, ('to_basis_set', ['.@'])
    )
    # add_mapping_annotation(basis_set.AtomCenteredBasisSet.basis_functions, OUT_KEY, '.x_lobster_basis.x_lobster_basis_species.2')


class BasisSetContainer(basis_set.BasisSetContainer):
    add_mapping_annotation(
        basis_set.AtomCenteredBasisSet.m_def, OUT_KEY, '.x_lobster_basis_species'
    )


class ModelMethod(model_method.ModelMethod):
    add_mapping_annotation(
        basis_set.BasisSetContainer.m_def, OUT_KEY, '.x_lobster_basis'
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.wall_start, OUT_KEY, ('to_unix_time', ['.datetime'])
    )
    add_mapping_annotation(general.Simulation.model_system, STRUCTURE_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_method, OUT_KEY, '.@')
    add_mapping_annotation(Outputs.m_def, ICOXPLIST_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, STRUCTURE_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, ICOXPLIST_KEY, '@')

m_package.__init_metainfo__()
