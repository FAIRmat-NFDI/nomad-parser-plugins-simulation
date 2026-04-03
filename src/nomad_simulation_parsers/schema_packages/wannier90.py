from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
    variables,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

WOUT_KEY = 'wannier_wout'
WIN_KEY = 'wannier_win'
BAND_KEY = 'wannier_band'
WHR_KEY = 'wannier_whr'
DOS_KEY = 'wannier_dos'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, WOUT_KEY, '.version')


class WannierSphericalSymmetryState(atoms_state.SphericalSymmetryState):
    """
    Spherical symmetry state customized for Wannier90 orbital projections.

    Maps Wannier90 input file orbital specifications to quantum numbers.
    """

    add_mapping_annotation(
        atoms_state.SphericalSymmetryState.l_quantum_number, WIN_KEY, '.l'
    )
    add_mapping_annotation(
        atoms_state.SphericalSymmetryState.ml_quantum_number, WIN_KEY, '.m'
    )


# class OrbitalsState(atoms_state.OrbitalsState):
#     atoms_state.OrbitalsState.l_quantum_symbol.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(win=Mapper(mapper='.l')))

#     atoms_state.OrbitalsState.ml_quantum_symbol.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(win=Mapper(mapper='.m')))


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.chemical_symbol.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(wout=Mapper(mapper='.@')))

    # model_system.AtomsState.orbitals_state.m_annotations.setdefault(
    #     MAPPING_ANNOTATION_KEY, {}
    # ).update(dict(win=Mapper(mapper=('get_orbitals_state', ['.projection[1]']))))


class AtomicCell(model_system.Representation):
    model_system.Representation.lattice_vectors.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            wout=Mapper(
                mapper=('get_lattice_vectors', ['lattice_vectors']), unit='angstrom'
            )
        )
    )


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors,
        WOUT_KEY,
        ('get_lattice_vectors', ['lattice_vectors']),
        unit='angstrom',
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        WOUT_KEY,
        ('get_pbc', ['lattice_vectors']),
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, WOUT_KEY, '.@')
    add_mapping_annotation(model_system.ModelSystem.representations, WIN_KEY, '.@')
    add_mapping_annotation(
        model_system.ModelSystem.sub_systems,
        WIN_KEY,
        ('get_projections', ['.projections']),
    )
    add_mapping_annotation(
        model_system.ModelSystem.branch_label,
        WIN_KEY,
        (
            'get_branch_label_indices',
            [
                '.projection[0]',
                'structure.positions',
                'structure.labels',
                'lattice_vectors',
            ],
        ),
        search='label',
        cache=True,
    )
    add_mapping_annotation(
        model_system.ModelSystem.particle_indices,
        WIN_KEY,
        (
            'get_branch_label_indices',
            [
                '.projection[0]',
                'structure.positions',
                'structure.labels',
                'lattice_vectors',
            ],
        ),
        search='indices',
        cache=True,
    )
    add_mapping_annotation(model_system.AtomsState.m_def, WOUT_KEY, '.labels')
    add_mapping_annotation(model_system.AtomsState.m_def, WIN_KEY, '.@')
    add_mapping_annotation(
        model_system.ModelSystem.positions, WOUT_KEY, '.positions', unit='angstrom'
    )


class KMesh(numerical_settings.KMesh):
    add_mapping_annotation(numerical_settings.KMesh.n_points, WOUT_KEY, '.n_points')
    add_mapping_annotation(numerical_settings.KMesh.grid, WOUT_KEY, '.grid')
    add_mapping_annotation(
        numerical_settings.KMesh.points, WOUT_KEY, ('get_kpoints', ['.k_points'])
    )


class KLinePath(numerical_settings.KLinePath):
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_names, WOUT_KEY, '.names'
    )
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_values, WOUT_KEY, '.values'
    )
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_names, BAND_KEY, '.names'
    )
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_values, BAND_KEY, '.values'
    )


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(numerical_settings.KSpace.k_mesh, WOUT_KEY, '.k_mesh')
    add_mapping_annotation(
        numerical_settings.KSpace.k_line_path,
        WOUT_KEY,
        ('get_k_line_path', ['.k_line_path']),
        cache=True,
    )


class ModelMethod(model_method.ModelMethod):
    add_mapping_annotation(numerical_settings.KSpace.m_def, WOUT_KEY, '.@')


class Wannier(model_method.Wannier):
    add_mapping_annotation(
        model_method.Wannier.is_maximally_localized,
        WOUT_KEY,
        ('is_maximally_localized', ['.Niter'], dict(default=0)),
    )
    add_mapping_annotation(
        model_method.Wannier.energy_window_outer, WOUT_KEY, '.energy_windows.outer'
    )
    add_mapping_annotation(
        model_method.Wannier.energy_window_inner, WOUT_KEY, '.energy_windows.inner'
    )
    # add_mapping_annotations(model_method.Wannier.n_orbitals, WOUT_KEY, '.Nwannier')


# TODO: check whether this section is k-dependent
class ElectronicBandStructure(outputs.ElectronicBandStructure):
    # properties.ElectronicBandStructure.n_bands.m_annotations.setdefault(
    #     MAPPING_ANNOTATION_KEY, {}
    # ).update(dict(wout=Mapper(mapper='.Nwannier')))

    add_mapping_annotation(
        outputs.ElectronicBandStructure.value, BAND_KEY, ('get_data', ['.data'])
    )
    add_mapping_annotation(outputs.ElectronicBandStructure.k_path, BAND_KEY, '.k_path')


class WignerSeitz(variables.WignerSeitz):
    add_mapping_annotation(variables.WignerSeitz.n_points, WHR_KEY, '.n_ws_points')
    add_mapping_annotation(variables.WignerSeitz.points, WHR_KEY, '.ws_points')


class HoppingMatrix(properties.HoppingMatrix):
    add_mapping_annotation(properties.HoppingMatrix.n_orbitals, WHR_KEY, '.n_orbitals')
    add_mapping_annotation(
        properties.HoppingMatrix.degeneracy_factors, WHR_KEY, '.degeneracy_factors'
    )

    # TODO shape mismatch
    # add_mapping_annotation(
    #     properties.HoppingMatrix.value, WHR_KEY, '.hoppings', unit='eV'
    # )

    add_mapping_annotation(variables.WignerSeitz.m_def, WHR_KEY, '.@')


class CrystalFieldSplitting(properties.CrystalFieldSplitting):
    add_mapping_annotation(
        properties.CrystalFieldSplitting.n_orbitals, WHR_KEY, 'n_orbitals'
    )
    add_mapping_annotation(
        properties.CrystalFieldSplitting.value, WHR_KEY, '.crystal_fields', unit='eV'
    )


class Energy2(variables.Energy2):
    add_mapping_annotation(variables.Energy2.points, DOS_KEY, '.energies', unit='eV')


class ElectronicDensityOfStates(properties.ElectronicDensityOfStates):
    add_mapping_annotation(
        properties.ElectronicDensityOfStates.value, DOS_KEY, '.value', unit='1/eV'
    )
    add_mapping_annotation(variables.Energy2.m_def, DOS_KEY, '.@')


class Outputs(outputs.Outputs):
    # Legacy parity: Wannier90 band structures come from `*band.dat`; avoid
    # placeholder sections from `.wout` metadata-only mappings.
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        BAND_KEY,
        ('get_band_structure', ['.data', '.k_path']),
    )
    add_mapping_annotation(
        outputs.Outputs.hopping_matrices,
        WHR_KEY,
        ('get_hoppings', ['.@'], dict(ws=True)),
    )
    add_mapping_annotation(
        outputs.Outputs.crystal_field_splittings, WHR_KEY, ('get_hoppings', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos, DOS_KEY, ('get_dos', ['.data'])
    )
    # TODO(legacy-parity): legacy Wannier90 parser did not populate explicit
    # electronic band-gap sections; keep unmapped for now.


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, WOUT_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_system, WIN_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_system, WOUT_KEY, '.structure')
    add_mapping_annotation(model_method.Wannier.m_def, WOUT_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, WOUT_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, BAND_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, WHR_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, DOS_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, WIN_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, WOUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, BAND_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, WHR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOS_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
