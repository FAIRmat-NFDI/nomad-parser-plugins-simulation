from nomad.metainfo import Quantity, SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

OUT_KEY = 'yambo_out'
NETCDF_KEY = 'yambo_netcdf'
SPECTRA_KEY = 'yambo_spectra'

m_package = SchemaPackage()


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.version')


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, NETCDF_KEY, '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        NETCDF_KEY,
        ('get_positions', []),
        unit='angstrom',
    )
    add_mapping_annotation(
        model_system.ModelSystem.lattice_vectors,
        NETCDF_KEY,
        ('get_lattice_vectors', []),
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def, NETCDF_KEY, ('get_labels', [])
    )


class KMesh(numerical_settings.KMesh):
    add_mapping_annotation(numerical_settings.KMesh.all_points, NETCDF_KEY, '.@')


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(
        numerical_settings.KSpace.k_mesh, NETCDF_KEY, ('get_kpoints', [])
    )


class ModelMethod(model_method.ModelMethod):
    add_mapping_annotation(numerical_settings.KSpace.m_def, NETCDF_KEY, '.@')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, NETCDF_KEY, '.energies')
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, OUT_KEY, '.energies')


class AbsorptionSpectra(outputs.AbsorptionSpectrum):
    add_mapping_annotation(
        outputs.AbsorptionSpectrum.value, SPECTRA_KEY, '.intensities'
    )
    add_mapping_annotation(
        outputs.AbsorptionSpectrum.energies, SPECTRA_KEY, '.excitation_energies'
    )
# EM, 6 Jul, 2026:    
    add_mapping_annotation(
        outputs.AbsorptionSpectrum.sp_type, SPECTRA_KEY, 'sp_type'
    )
    add_mapping_annotation(
        outputs.AbsorptionSpectrum.n_energies, SPECTRA_KEY, 'n_energies'
    )
# end EM


class Outputs(outputs.Outputs):
    # TODO add description
#    sp_type = Quantity(type=str)   # EM: commented,  6 Jul, 2026

#    add_mapping_annotation(sp_type, SPECTRA_KEY, 'sp_type')  # EM: changed sp_type to SPECTRA_KEY, Jul 1st, 2026; commented whole line:  6 Jul, 2026

    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        NETCDF_KEY,
        ('get_eigenvalues', []),
        cache=True,
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues, OUT_KEY, '.eigenvalues'
    )
    # TODO is this the correct def to use for spectra data
    add_mapping_annotation(
        AbsorptionSpectrum.m_def, SPECTRA_KEY, ('get_spectra', ['data'])  # EM: outputs.AbsorptionSpectrum  -->  AbsorptionSpectrum
    )


class Simulation(general.Simulation):
    add_mapping_annotation(
        general.Simulation.wall_start, OUT_KEY, ('get_wallstart', ['.date_start'])
    )
    add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_system, NETCDF_KEY, '.@')
    add_mapping_annotation(general.Simulation.model_method, NETCDF_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, NETCDF_KEY, '.@')
    add_mapping_annotation(
        Outputs.m_def,
        OUT_KEY,
        (
            'get_outputs',
            [
                '.core_variables_setup.energies_occupations',
                '.modules',
                '.transferred_momenta',
            ],
        ),
    )
    add_mapping_annotation(Outputs.m_def, SPECTRA_KEY, '.@')

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, NETCDF_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, SPECTRA_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    # Metainfo initialization errors are intentionally ignored here to avoid
    # failing on import in environments where the NOMAD metainfo registry or
    # plugin infrastructure is not fully configured.
    pass
