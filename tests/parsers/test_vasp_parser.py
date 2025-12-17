from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

LOGGER = get_logger(__name__)


def test_vasprun():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, LOGGER)

    # Check numerical_settings structure
    simulation = archive.data
    assert simulation is not None

    LOGGER.info('\n=== vasprun.xml structure ===')
    if hasattr(simulation, 'model_method') and simulation.model_method:
        for method in simulation.model_method:
            LOGGER.info(f'ModelMethod type: {type(method).__name__}')
            if hasattr(method, 'numerical_settings') and method.numerical_settings:
                LOGGER.info(
                    f'Found {len(method.numerical_settings)} numerical_settings'
                )
                for ns in method.numerical_settings:
                    LOGGER.info(f'  - {type(ns).__name__}')
                    if hasattr(ns, 'pseudopotential') and ns.pseudopotential:
                        LOGGER.info(
                            f'    Contains {len(ns.pseudopotential)} pseudopotentials'
                        )
                        for pp in ns.pseudopotential:
                            LOGGER.info(f'      * {pp.name}')


def test_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)

    # Verify numerical_settings are present
    simulation = archive.data
    assert simulation is not None, 'No simulation data in archive'
    assert simulation.model_method, 'No model_method in simulation'

    method = simulation.model_method[0]
    assert hasattr(method, 'numerical_settings'), (
        'ModelMethod missing numerical_settings attribute'
    )
    assert method.numerical_settings is not None, 'numerical_settings is None'
    assert len(method.numerical_settings) > 0, 'numerical_settings is empty'

    LOGGER.info(f'Found {len(method.numerical_settings)} numerical_settings')

    # Check for Pseudopotentials
    pseudopotentials = [
        ns for ns in method.numerical_settings if type(ns).__name__ == 'Pseudopotential'
    ]
    assert len(pseudopotentials) > 0, 'No Pseudopotential objects in numerical_settings'

    LOGGER.info(f'Found {len(pseudopotentials)} Pseudopotentials:')
    for pp in pseudopotentials:
        LOGGER.info(f'  - {pp.name}: {pp.n_valence_electrons} valence electrons')
        # Verify OUTCAR-specific fields are populated
        if pp.sha256:
            LOGGER.info(f'    SHA256: {pp.sha256[:16]}...')
        if pp.l_max is not None:
            LOGGER.info(f'    l_max: {pp.l_max}')
