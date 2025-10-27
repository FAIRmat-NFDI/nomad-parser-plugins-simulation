from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

INFO_KEY = 'info'
INPUT_XML_KEY = 'input_xml'
EIGVAL_KEY = 'eigval'
BANDSTRUCTURE_XML_KEY = 'bandstructure_xml'
DOS_XML_KEY = 'dos_xml'

# simulation
add_mapping_annotation(general.Simulation.m_def, INFO_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, INPUT_XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EIGVAL_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, BANDSTRUCTURE_XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOS_XML_KEY, '@')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, INFO_KEY, '.@')
    # DFT method
    add_mapping_annotation(
        model_method.DFT.m_def, INFO_KEY, '.initialization.xc_functional'
    )
    add_mapping_annotation(model_method.DFT.m_def, INPUT_XML_KEY, '.input.groundstate')
    add_mapping_annotation(model_method.DFT.m_def, BANDSTRUCTURE_XML_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.model_system,
        INFO_KEY,
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotation(
        general.Simulation.outputs, INFO_KEY, ('get_configurations', ['.@'])
    )
    add_mapping_annotation(general.Simulation.outputs, EIGVAL_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, BANDSTRUCTURE_XML_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, DOS_XML_KEY, '.@')


class Program(general.Program):
    add_mapping_annotation(general.Program.version, INFO_KEY, '.program_version')


class ModelMethod(model_method.ModelMethod):
    add_mapping_annotation(
        numerical_settings.KSpace.m_def, BANDSTRUCTURE_XML_KEY, '.@'
    )


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(
        numerical_settings.KSpace.k_line_path, BANDSTRUCTURE_XML_KEY, '.@'
    )


class KLinePath(numerical_settings.KLinePath):
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_names,
        BANDSTRUCTURE_XML_KEY,
        r'bandstructure.vertex[*]."@label"',
    )
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_values,
        BANDSTRUCTURE_XML_KEY,
        ('reshape_coords', [r'bandstructure.vertex[*]."@coord"']),
    )


class DFT(model_method.DFT):
    add_mapping_annotation(
        model_method.DFT.xc_functionals, INFO_KEY, ('get_xc_functionals', ['.type'])
    )
    add_mapping_annotation(
        model_method.DFT.xc_functionals,
        INPUT_XML_KEY,
        ('get_xc_functionals', ['.libxc']),
    )


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(model_method.XCFunctional.libxc_name, INFO_KEY, '.libxc')
    add_mapping_annotation(
        model_method.XCFunctional.libxc_name, INPUT_XML_KEY, '.libxc'
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.AtomicCell.m_def, INFO_KEY, '.@')
    add_mapping_annotation(model_system.ModelSystem.positions, INFO_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, INFO_KEY, '.atoms')


class AtomicCell(model_system.AtomicCell):
    add_mapping_annotation(
        model_system.AtomicCell.lattice_vectors, INFO_KEY, '.lattice_vectors'
    )


class AtomsState(atoms_state.AtomsState):
    add_mapping_annotation(atoms_state.AtomsState.chemical_symbol, INFO_KEY, '.symbol')


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, INFO_KEY, '.@')
    add_mapping_annotation(
        outputs.Outputs.total_forces, INFO_KEY, ('get_forces', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues, EIGVAL_KEY, ('get_eigenvalues', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        BANDSTRUCTURE_XML_KEY,
        ('get_bandstructures', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos, DOS_XML_KEY, 'dos.totaldos.diagram'
    )


class TotalEnergy(properties.TotalEnergy):
    add_mapping_annotation(
        properties.TotalEnergy.value, INFO_KEY, '.final.energy_total || energy_total'
    )


class TotalForce(properties.forces.TotalForce):
    add_mapping_annotation(properties.forces.TotalForce.value, INFO_KEY, '.forces')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.n_bands, EIGVAL_KEY, '.n_states'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, EIGVAL_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, EIGVAL_KEY, '.occupancies'
    )


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    add_mapping_annotation(
        outputs.ElectronicBandStructure.n_bands, BANDSTRUCTURE_XML_KEY, '.n_states'
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.value, BANDSTRUCTURE_XML_KEY, '.energies'
    )


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    ###### TODO read unit from axis
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value,
        DOS_XML_KEY,
        ('to_float', [r'.point[*]."@dos"']),
        unit='1/hartree',
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.projected_dos,
        DOS_XML_KEY,
        'dos.partialdos.diagram',
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
