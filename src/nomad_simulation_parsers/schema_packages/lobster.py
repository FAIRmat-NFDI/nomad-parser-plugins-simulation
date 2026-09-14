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
    variables,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'lobster_out'
STRUCTURE_KEY = 'lobster_structure'
ICOXPLIST_KEY = 'lobster_icoxplist'
COXPCAR_KEY = 'lobster_coxpcar'
CHARGE_KEY = 'lobster_charge'
DOSCAR_KEY = 'lobster_doscar'


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
    add_mapping_annotation(
        x_lobster_number_of_cohp_pairs,
        ICOXPLIST_KEY,
        ('get_atom_pair_indices', ['.@']),
        search='length(@)',
    )

    x_lobster_cohp_atom1_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Species and indices of the first atom for which is the specific COHP/iCOHP
        calculated
        """,
    )
    add_mapping_annotation(
        x_lobster_cohp_atom1_labels,
        ICOXPLIST_KEY,
        ('get_atom_labels', ['.@'], dict(atom=0)),
    )

    x_lobster_cohp_atom2_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Species and indices of the second atom for which is the specific COHP/iCOHP
        calculated
        """,
    )
    add_mapping_annotation(
        x_lobster_cohp_atom2_labels,
        ICOXPLIST_KEY,
        ('get_atom_labels', ['.@'], dict(atom=1)),
    )

    x_lobster_cohp_distances = Quantity(
        type=np.float64,
        unit='meter',
        shape=['x_lobster_number_of_cohp_pairs'],
        description="""
        Distance between atoms of the pair for which is the specific COHP/iCOHP
        calculated.
        """,
    )
    add_mapping_annotation(
        x_lobster_cohp_distances,
        ICOXPLIST_KEY,
        ('get_atom_distances', ['.@']),
        unit='angstrom',
    )

    x_lobster_cohp_translations = Quantity(
        type=np.int32,
        shape=['x_lobster_number_of_cohp_pairs', 3],
        description="""
        Vector connecting the unit-cell of the first atom with the one of the second
        atom

        This is only used with LOBSTER versions 3.0.0 and above, older versions use
        x_lobster_cohp_number_of_bonds instead.
        """,
    )
    add_mapping_annotation(
        x_lobster_cohp_translations, ICOXPLIST_KEY, ('get_translations', ['.@'])
    )

    x_lobster_integrated_cohp_at_fermi_level = Quantity(
        type=np.float64,
        unit='joule',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_pairs'],
        description="""
        Calculated iCOHP values integrated up to the Fermi level.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_cohp_at_fermi_level,
        ICOXPLIST_KEY,
        ('get_integrated_coxp_at_fermi_level', ['.@']),
        unit='eV',
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
    add_mapping_annotation(x_lobster_cohp_energies, COXPCAR_KEY, '.energy', unit='eV')

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
    add_mapping_annotation(
        x_lobster_cohp_values,
        COXPCAR_KEY,
        ('get_atom_value', ['.pair_coxp', '.coxp_pairs']),
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
    add_mapping_annotation(
        x_lobster_integrated_cohp_values,
        COXPCAR_KEY,
        ('get_atom_value', ['.pair_icoxp', '.coxp_pairs']),
        unit='eV',
    )

    x_lobster_average_cohp_values = Quantity(
        type=np.float64,
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_values'],
        description="""
        Calculated COHP values averaged over all pairs.
        """,
    )
    add_mapping_annotation(x_lobster_average_cohp_values, COXPCAR_KEY, '.total_coxp')

    x_lobster_average_integrated_cohp_values = Quantity(
        type=np.float64,
        unit='joule',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_values'],
        description="""
        Calculated iCOHP values averaged over all pairs.
        """,
    )
    add_mapping_annotation(
        x_lobster_average_integrated_cohp_values, COXPCAR_KEY, '.total_icoxp', unit='eV'
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
    add_mapping_annotation(
        x_lobster_cohp_number_of_bonds, ICOXPLIST_KEY, ('get_bonds', ['.@'])
    )

    x_lobster_cohp_orbital_per_label = SubSection(
        sub_section=SectionProxy('x_lobster_section_cohp_orbital_label'), repeats=True
    )
    add_mapping_annotation(x_lobster_cohp_orbital_per_label, COXPCAR_KEY, '.bond_pairs')
    add_mapping_annotation(
        x_lobster_cohp_orbital_per_label,
        ICOXPLIST_KEY,
        ('get_atom_pair_indices', ['.@']),
    )


class x_lobster_section_cohp_orbital_label(MSection):
    """
    Orbital-resolved COHP/iCOHP grouped per bond label.
    """

    m_def = Section(validate=False)

    x_lobster_pair_label = Quantity(
        type=str,
        description="""
        Bond label as printed in COHPCAR.
        """,
    )
    add_mapping_annotation(
        x_lobster_pair_label, COXPCAR_KEY, ('get_bond_label', ['.@'])
    )
    add_mapping_annotation(
        x_lobster_pair_label, ICOXPLIST_KEY, ('get_bond_label', ['.@', 'COHPCAR'])
    )

    x_lobster_orbital_pairs = SubSection(
        sub_section=SectionProxy('x_lobster_section_cohp_orbital_pair'), repeats=True
    )
    add_mapping_annotation(
        x_lobster_orbital_pairs,
        COXPCAR_KEY,
        ('get_orbital_pairs', ['.@', 'COHPCAR'], dict(name='cohp')),
    )


class x_lobster_section_cohp_orbital_pair(MSection):
    """
    Individual orbital pair data for COHP/iCOHP.
    """

    m_def = Section(validate=False)

    x_lobster_atom1_orbital = Quantity(
        type=str,
        description="""
        Orbital name for atom1 (e.g., 'Na1_3s').
        """,
    )
    add_mapping_annotation(
        x_lobster_atom1_orbital,
        COXPCAR_KEY,
        ('get_atom_label', ['.@', 'COHPCAR'], dict(atom=0)),
    )
    # add_mapping_annotation(x_lobster_atom1_orbital, ICOXPLIST_KEY, '.atomMU')

    x_lobster_atom2_orbital = Quantity(
        type=str,
        description="""
        Orbital name for atom2 (e.g., 'Cl2_3p').
        """,
    )
    add_mapping_annotation(
        x_lobster_atom2_orbital,
        COXPCAR_KEY,
        ('get_atom_label', ['.@', 'COHPCAR'], dict(atom=1)),
    )
    # add_mapping_annotation(x_lobster_atom2_orbital, ICOXPLIST_KEY, '.atomNU')

    x_lobster_cohp_orbital_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_values'],
        description="""
        Calculated COHP values for this specific orbital pair.
        """,
    )
    add_mapping_annotation(
        x_lobster_cohp_orbital_values,
        COXPCAR_KEY,
        ('get_orbital_value', ['.@', 'COHPCAR'], dict(type='pair_coxp')),
    )

    x_lobster_integrated_cohp_orbital_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cohp_values'],
        description="""
        Calculated iCOHP values for this specific orbital pair.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_cohp_orbital_values,
        COXPCAR_KEY,
        ('get_orbital_value', ['.@', 'COHPCAR'], dict(type='pair_icoxp', name='cohp')),
    )

    x_lobster_integrated_orbital_cohp_at_fermi_level = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels'],
        description="""
        Calculated orbital iCOHP values integrated up to the Fermi level for this
        specific orbital pair.
        """,
    )


class x_lobster_section_coop_orbital_label(MSection):
    """
    Orbital-resolved COOP/iCOOP grouped per bond label.
    """

    m_def = Section(validate=False)

    x_lobster_pair_label = Quantity(
        type=str,
        description="""
        Bond label as printed in COOPCAR.
        """,
    )

    add_mapping_annotation(
        x_lobster_pair_label, COXPCAR_KEY, ('get_bond_label', ['.@'])
    )

    x_lobster_orbital_pairs = SubSection(
        sub_section=SectionProxy('x_lobster_section_coop_orbital_pair'), repeats=True
    )
    add_mapping_annotation(
        x_lobster_orbital_pairs, COXPCAR_KEY, ('get_orbital_pairs', ['.@', 'COOPCAR'])
    )


class x_lobster_section_coop_orbital_pair(MSection):
    """
    Individual orbital pair data for COOP/iCOOP.
    """

    m_def = Section(validate=False)

    x_lobster_atom1_orbital = Quantity(
        type=str,
        description="""
        Orbital name for atom1 (e.g., 'Na1_3s').
        """,
    )
    add_mapping_annotation(
        x_lobster_atom1_orbital,
        COXPCAR_KEY,
        ('get_atom_label', ['.@', 'COHPCAR'], dict(atom=0)),
    )
    # add_mapping_annotation(x_lobster_atom1_orbital, ICOXPLIST_KEY, '.atomMU')

    x_lobster_atom2_orbital = Quantity(
        type=str,
        description="""
        Orbital name for atom2 (e.g., 'Cl2_3p').
        """,
    )
    add_mapping_annotation(
        x_lobster_atom2_orbital,
        COXPCAR_KEY,
        ('get_atom_label', ['.@', 'COHPCAR'], dict(atom=1)),
    )
    # add_mapping_annotation(x_lobster_atom2_orbital, ICOXPLIST_KEY, '.atomNU')

    x_lobster_coop_orbital_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_coop_values'],
        description="""
        Calculated COOP values for this specific orbital pair.
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_orbital_values,
        COXPCAR_KEY,
        ('get_orbital_value', ['.@', 'COOPCAR'], dict(type='pair_coxp')),
    )

    x_lobster_integrated_coop_orbital_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_coop_values'],
        description="""
        Calculated iCOOP values for this specific orbital pair.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_coop_orbital_values,
        COXPCAR_KEY,
        ('get_orbital_value', ['.@', 'COOPCAR'], dict(type='pair_icoxp')),
    )

    x_lobster_integrated_orbital_coop_at_fermi_level = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels'],
        description="""
        Calculated orbital iCOOP values integrated up to the Fermi level for this
        specific orbital pair.
        """,
    )


class x_lobster_section_coop(MSection):
    """
    This is a section containing the crystal orbital hamilton population (COOP)
    and integrated coop (iCOOP) values.
    """

    m_def = Section(validate=False)

    x_lobster_number_of_coop_pairs = Quantity(
        type=int,
        description="""
        Number of atom pairs for which are the COOPs and iCOOPs calculated.
        """,
    )
    add_mapping_annotation(
        x_lobster_number_of_coop_pairs,
        ICOXPLIST_KEY,
        ('get_atom_pair_indices', ['.@']),
        search='length(@)',
    )

    x_lobster_coop_atom1_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_coop_pairs'],
        description="""
        Species and indices of the first atom for which is the specific COOP/iCOOP
        calculated
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_atom1_labels,
        ICOXPLIST_KEY,
        ('get_atom_labels', ['.@'], dict(atom=0)),
    )

    x_lobster_coop_atom2_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_coop_pairs'],
        description="""
        Species and indices of the second atom for which is the specific COOP/iCOOP
        calculated
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_atom2_labels,
        ICOXPLIST_KEY,
        ('get_atom_labels', ['.@'], dict(atom=1)),
    )

    x_lobster_coop_distances = Quantity(
        type=np.float64,
        unit='meter',
        shape=['x_lobster_number_of_coop_pairs'],
        description="""
        Distance between atoms of the pair for which is the specific COOP/iCOOP
        calculated.
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_distances,
        ICOXPLIST_KEY,
        ('get_atom_distances', ['.@']),
        unit='angstrom',
    )

    x_lobster_coop_translations = Quantity(
        type=np.int32,
        shape=['x_lobster_number_of_coop_pairs', 3],
        description="""
        Vector connecting the unit-cell of the first atom with the one of the second
        atom

        This is only used with LOBSTER versions 3.0.0 and above, older versions use
        x_lobster_coop_number_of_bonds instead.
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_translations, ICOXPLIST_KEY, ('get_translations', ['.@'])
    )

    x_lobster_integrated_coop_at_fermi_level = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_coop_pairs'],
        description="""
        Calculated iCOOP values integrated up to the Fermi level.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_coop_at_fermi_level,
        ICOXPLIST_KEY,
        ('get_integrated_coxp_at_fermi_level', ['.@']),
    )

    x_lobster_number_of_coop_values = Quantity(
        type=int,
        description="""
        Number of energy values for the COOP and iCOOP.
        """,
    )

    x_lobster_coop_energies = Quantity(
        type=np.float64,
        unit='joule',
        shape=['x_lobster_number_of_coop_values'],
        description="""
        Array containing the set of discrete energy values for COOP and iCOOP.
        """,
    )
    add_mapping_annotation(x_lobster_coop_energies, COXPCAR_KEY, '.energy', unit='eV')

    x_lobster_coop_values = Quantity(
        type=np.float64,
        shape=[
            'x_lobster_number_of_coop_pairs',
            'number_of_spin_channels',
            'x_lobster_number_of_coop_values',
        ],
        description="""
        Calculated COOP values.
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_values,
        COXPCAR_KEY,
        ('get_atom_value', ['.pair_coxp', '.coxp_pairs']),
    )

    x_lobster_integrated_coop_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=[
            'x_lobster_number_of_coop_pairs',
            'number_of_spin_channels',
            'x_lobster_number_of_coop_values',
        ],
        description="""
        Calculated iCOOP values.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_coop_values,
        COXPCAR_KEY,
        ('get_atom_value', ['.pair_icoxp', '.coxp_pairs']),
    )

    x_lobster_average_coop_values = Quantity(
        type=np.float64,
        shape=['number_of_spin_channels', 'x_lobster_number_of_coop_values'],
        description="""
        Calculated COOP values averaged over all pairs.
        """,
    )
    add_mapping_annotation(x_lobster_average_coop_values, COXPCAR_KEY, '.total_coxp')

    x_lobster_average_integrated_coop_values = Quantity(
        type=np.float32,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_coop_values'],
        description="""
        Calculated iCOOP values averaged over all pairs.
        """,
    )
    add_mapping_annotation(
        x_lobster_average_integrated_coop_values, COXPCAR_KEY, '.total_icoxp'
    )

    x_lobster_coop_number_of_bonds = Quantity(
        type=int,
        shape=['x_lobster_number_of_coop_pairs'],
        description="""
        Number of bonds between first atom and the second atom (including
        the periodic images).

        This is only used in older LOBSTER versions, new versions print one line
        for every neighbor, so a pair which had x_lobster_icoop_number_of_bonds = 4
        in the old version would actually show as 4 lines in the ICOOPLIST or 4 columns
        in the COOPCAR in the new format.
        """,
    )
    add_mapping_annotation(
        x_lobster_coop_number_of_bonds, ICOXPLIST_KEY, ('get_bonds', ['.@'])
    )

    x_lobster_coop_orbital_per_label = SubSection(
        sub_section=SectionProxy('x_lobster_section_coop_orbital_label'), repeats=True
    )
    add_mapping_annotation(x_lobster_coop_orbital_per_label, COXPCAR_KEY, '.bond_pairs')
    add_mapping_annotation(
        x_lobster_coop_orbital_per_label,
        ICOXPLIST_KEY,
        ('get_atom_pair_indices', ['.@']),
    )


class x_lobster_section_cobi(MSection):
    """
    This is a section containing the crystal orbital bond index (COBI)
    and integrated cobi (iCOBI) values.
    """

    m_def = Section(validate=False)

    x_lobster_number_of_cobi_pairs = Quantity(
        type=int,
        description="""
        Number of atom pairs for which are the COBIs and iCOBIs calculated.
        """,
    )
    add_mapping_annotation(
        x_lobster_number_of_cobi_pairs,
        ICOXPLIST_KEY,
        ('get_atom_pair_indices', ['.@']),
        search='length(@)',
    )

    x_lobster_cobi_atom1_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_cobi_pairs'],
        description="""
        Species and indices of the first atom for which is the specific COBI/iCOBI
        calculated
        """,
    )
    add_mapping_annotation(
        x_lobster_cobi_atom1_labels,
        ICOXPLIST_KEY,
        ('get_atom_labels', ['.@'], dict(atom=0)),
    )

    x_lobster_cobi_atom2_labels = Quantity(
        type=str,
        shape=['x_lobster_number_of_cobi_pairs'],
        description="""
        Species and indices of the second atom for which is the specific COBI/iCOBI
        calculated
        """,
    )
    add_mapping_annotation(
        x_lobster_cobi_atom2_labels,
        ICOXPLIST_KEY,
        ('get_atom_labels', ['.@'], dict(atom=1)),
    )

    x_lobster_cobi_distances = Quantity(
        type=np.float64,
        unit='meter',
        shape=['x_lobster_number_of_cobi_pairs'],
        description="""
        Distance between atoms of the pair for which is the specific COBI/iCOBI
        calculated.
        """,
    )
    add_mapping_annotation(
        x_lobster_cobi_distances,
        ICOXPLIST_KEY,
        ('get_atom_distances', ['.@']),
        unit='angstrom',
    )

    x_lobster_cobi_translations = Quantity(
        type=np.int32,
        shape=['x_lobster_number_of_cobi_pairs', 3],
        description="""
        Vector connecting the unit-cell of the first atom with the one of the second
        atom

        This is only used with LOBSTER versions 4.1.0 and above.
        """,
    )
    add_mapping_annotation(
        x_lobster_cobi_translations, ICOXPLIST_KEY, ('get_translations', ['.@'])
    )

    x_lobster_integrated_cobi_at_fermi_level = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cobi_pairs'],
        description="""
        Calculated iCOBI values integrated up to the Fermi level.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_cobi_at_fermi_level,
        ICOXPLIST_KEY,
        ('get_integrated_coxp_at_fermi_level', ['.@']),
        unit='eV',
    )

    x_lobster_number_of_cobi_values = Quantity(
        type=int,
        description="""
        Number of energy values for the COBI and iCOBI.
        """,
    )

    x_lobster_cobi_energies = Quantity(
        type=np.float64,
        unit='joule',
        shape=['x_lobster_number_of_cobi_values'],
        description="""
        Array containing the set of discrete energy values for COBI and iCOBI.
        """,
    )
    add_mapping_annotation(x_lobster_cobi_energies, COXPCAR_KEY, '.energy', unit='eV')

    x_lobster_cobi_values = Quantity(
        type=np.float64,
        shape=[
            'x_lobster_number_of_cobi_pairs',
            'number_of_spin_channels',
            'x_lobster_number_of_cobi_values',
        ],
        description="""
        Calculated COBI values.
        """,
    )
    add_mapping_annotation(
        x_lobster_cobi_values,
        COXPCAR_KEY,
        ('get_atom_value', ['.pair_coxp', '.coxp_pairs']),
    )

    x_lobster_integrated_cobi_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=[
            'x_lobster_number_of_cobi_pairs',
            'number_of_spin_channels',
            'x_lobster_number_of_cobi_values',
        ],
        description="""
        Calculated iCOBI values.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_cobi_values,
        COXPCAR_KEY,
        ('get_atom_value', ['.pair_icoxp', '.coxp_pairs']),
    )

    x_lobster_average_cobi_values = Quantity(
        type=np.float64,
        shape=['number_of_spin_channels', 'x_lobster_number_of_cobi_values'],
        unit='dimensionless',
        description="""
                Calculated COBI values averaged over all pairs.
                """,
    )
    add_mapping_annotation(x_lobster_average_cobi_values, COXPCAR_KEY, '.total_coxp')

    x_lobster_average_integrated_cobi_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cobi_values'],
        description="""
        Calculated iCOBI values averaged over all pairs.
        """,
    )
    add_mapping_annotation(
        x_lobster_average_integrated_cobi_values, COXPCAR_KEY, '.total_icoxp'
    )

    x_lobster_cobi_orbital_per_label = SubSection(
        sub_section=SectionProxy('x_lobster_section_cobi_orbital_label'), repeats=True
    )
    add_mapping_annotation(x_lobster_cobi_orbital_per_label, COXPCAR_KEY, '.bond_pairs')
    add_mapping_annotation(
        x_lobster_cobi_orbital_per_label,
        ICOXPLIST_KEY,
        ('get_atom_pair_indices', ['.@']),
    )


class x_lobster_section_cobi_orbital_label(MSection):
    """
    Orbital-resolved COBI/iCOBI grouped per bond label.
    """

    m_def = Section(validate=False)

    x_lobster_pair_label = Quantity(
        type=str,
        description="""
        Bond label as printed in COBICAR.
        """,
    )
    add_mapping_annotation(
        x_lobster_pair_label, COXPCAR_KEY, ('get_bond_label', ['.@'])
    )

    x_lobster_orbital_pairs = SubSection(
        sub_section=SectionProxy('x_lobster_section_cobi_orbital_pair'), repeats=True
    )
    add_mapping_annotation(
        x_lobster_orbital_pairs,
        COXPCAR_KEY,
        ('get_orbital_pairs', ['.@', 'COBICAR'], dict(name='COBI')),
    )


class x_lobster_section_cobi_orbital_pair(MSection):
    """
    Individual orbital pair data for COBI/iCOBI.
    """

    m_def = Section(validate=False)

    x_lobster_atom1_orbital = Quantity(
        type=str,
        description="""
        Orbital name for atom1 (e.g., 'Na1_3s').
        """,
    )
    add_mapping_annotation(
        x_lobster_atom1_orbital,
        COXPCAR_KEY,
        ('get_atom_label', ['.@', 'COHPCAR'], dict(atom=0)),
    )
    # add_mapping_annotation(x_lobster_atom1_orbital, ICOXPLIST_KEY, '.atomMU')

    x_lobster_atom2_orbital = Quantity(
        type=str,
        description="""
        Orbital name for atom2 (e.g., 'Cl2_3p').
        """,
    )
    add_mapping_annotation(
        x_lobster_atom2_orbital,
        COXPCAR_KEY,
        ('get_atom_label', ['.@', 'COHPCAR'], dict(atom=1)),
    )
    # add_mapping_annotation(x_lobster_atom2_orbital, ICOXPLIST_KEY, '.atomNU')

    x_lobster_cobi_orbital_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cobi_values'],
        description="""
        Calculated COBI values for this specific orbital pair.
        """,
    )
    add_mapping_annotation(
        x_lobster_cobi_orbital_values,
        COXPCAR_KEY,
        ('get_orbital_value', ['.@', 'COBICAR'], dict(type='pair_coxp', name='cobi')),
    )

    x_lobster_integrated_cobi_orbital_values = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels', 'x_lobster_number_of_cobi_values'],
        description="""
        Calculated iCOBI values for this specific orbital pair.
        """,
    )
    add_mapping_annotation(
        x_lobster_integrated_cobi_orbital_values,
        COXPCAR_KEY,
        ('get_orbital_value', ['.@', 'COBICAR'], dict(type='pair_icoxp', name='cobi')),
    )

    x_lobster_integrated_orbital_cobi_at_fermi_level = Quantity(
        type=np.float64,
        unit='dimensionless',
        shape=['number_of_spin_channels'],
        description="""
        Calculated orbital iCOBI values integrated up to the Fermi level for this
        specific orbital pair.
        """,
    )


# TODO migrate to general outputs section
class x_lobster_section_charge(MSection):
    m_def = Section(validate=False)

    label = Quantity(
        type=str,
    )
    add_mapping_annotation(label, CHARGE_KEY, '.symbol')

    type = Quantity(
        type=str,
    )
    add_mapping_annotation(type, CHARGE_KEY, '.kind')

    value = Quantity(type=np.float64, unit='coulomb')
    add_mapping_annotation(value, CHARGE_KEY, '.value', unit='elementary_charge')

    contributions = SubSection(
        sub_section=SectionProxy('x_lobster_section_charge'), repeats=True
    )
    add_mapping_annotation(contributions, CHARGE_KEY, '.contributions')


class Energy2(variables.Energy2):
    add_mapping_annotation(
        variables.Energy2.points, DOSCAR_KEY, 'DOSCAR.energies', unit='eV'
    )


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.energies, DOSCAR_KEY, 'DOSCAR', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value, DOSCAR_KEY, '.dos', unit='1/eV'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.projected_dos, DOSCAR_KEY, '.projected'
    )


class Outputs(outputs.Outputs):
    x_lobster_abs_total_spilling = Quantity(
        type=float,
        shape=['number_of_spin_channels'],
        description="""
        Absolute total spilling (in all levels)
        when projecting from the original wave functions into the local basis.
        """,
    )
    add_mapping_annotation(
        x_lobster_abs_total_spilling,
        OUT_KEY,
        ('get_spilling', ['.spilling'], dict(type='total')),
    )

    x_lobster_abs_charge_spilling = Quantity(
        type=float,
        shape=['number_of_spin_channels'],
        description="""
        Absolute total spilling of density (in occupied levels)
        when projecting from the original wave functions into the local basis.
        """,
    )
    add_mapping_annotation(
        x_lobster_abs_charge_spilling,
        OUT_KEY,
        ('get_spilling', ['.spilling'], dict(type='charge')),
    )

    x_lobster_section_cohp = SubSection(
        sub_section=SectionProxy('x_lobster_section_cohp')
    )
    add_mapping_annotation(x_lobster_section_cohp, ICOXPLIST_KEY, 'ICOHPLIST')
    add_mapping_annotation(x_lobster_section_cohp, COXPCAR_KEY, 'COHPCAR')

    x_lobster_section_coop = SubSection(
        sub_section=SectionProxy('x_lobster_section_coop')
    )
    add_mapping_annotation(x_lobster_section_coop, ICOXPLIST_KEY, 'ICOOPLIST')
    add_mapping_annotation(x_lobster_section_coop, COXPCAR_KEY, 'COOPCAR')

    x_lobster_section_cobi = SubSection(
        sub_section=SectionProxy('x_lobster_section_cobi')
    )
    add_mapping_annotation(x_lobster_section_cobi, ICOXPLIST_KEY, 'ICOBILIST')
    add_mapping_annotation(x_lobster_section_cobi, COXPCAR_KEY, 'COBICAR')

    x_lobster_section_charges = SubSection(
        sub_section=SectionProxy('x_lobster_section_charge'), repeats=True
    )
    add_mapping_annotation(
        x_lobster_section_charges, CHARGE_KEY, ('get_charges', ['.CHARGE'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos,
        DOSCAR_KEY,
        ('get_dos', ['.total_dos', '.projected_dos']),
    )


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(
        model_system.AtomsState.atomic_number, STRUCTURE_KEY, '.number'
    )
    add_mapping_annotation(model_system.AtomsState.atomic_number, DOSCAR_KEY, '.@')
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
    add_mapping_annotation(model_system.AtomsState.m_def, DOSCAR_KEY, '.atomic_numbers')
    add_mapping_annotation(
        general.ModelSystem.lattice_vectors, STRUCTURE_KEY, '.cell', unit='angstrom'
    )
    add_mapping_annotation(
        general.ModelSystem.periodic_boundary_conditions, STRUCTURE_KEY, '.pbc'
    )
    add_mapping_annotation(
        general.ModelSystem.periodic_boundary_conditions, DOSCAR_KEY, '.pbc'
    )


class AtomCenteredBasisSet(basis_set.AtomCenteredBasisSet):
    add_mapping_annotation(
        basis_set.AtomCenteredBasisSet.basis_set, OUT_KEY, ('to_basis_set', ['.@'])
    )


class BasisSetContainer(basis_set.BasisSetContainer):
    add_mapping_annotation(
        basis_set.AtomCenteredBasisSet.m_def, OUT_KEY, '.x_lobster_basis_species'
    )


class ModelMethod(model_method.ModelMethod):
    x_lobster_code = Quantity(
        type=str,
        description="""
        Used PAW program
        """,
    )
    add_mapping_annotation(x_lobster_code, OUT_KEY, '.x_lobster_code')

    add_mapping_annotation(
        basis_set.BasisSetContainer.m_def, OUT_KEY, '.x_lobster_basis'
    )


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.wall_start, OUT_KEY, ('to_unix_time', ['.datetime'])
    )
    add_mapping_annotation(general.Simulation.model_system, STRUCTURE_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_system, DOSCAR_KEY, '.DOSCAR')
    add_mapping_annotation(ModelMethod.m_def, OUT_KEY, '.@')
    add_mapping_annotation(Outputs.m_def, OUT_KEY, '.@')
    add_mapping_annotation(Outputs.m_def, ICOXPLIST_KEY, '.@')
    add_mapping_annotation(Outputs.m_def, COXPCAR_KEY, '.@')
    add_mapping_annotation(Outputs.m_def, CHARGE_KEY, '.@')
    add_mapping_annotation(Outputs.m_def, DOSCAR_KEY, '.DOSCAR')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, STRUCTURE_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, ICOXPLIST_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, COXPCAR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, CHARGE_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOSCAR_KEY, '@')

m_package.__init_metainfo__()
