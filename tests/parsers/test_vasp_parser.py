from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

LOGGER = get_logger(__name__)


def test_vasprun():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, LOGGER)

    simulation = archive.data
    assert simulation is not None


def test_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)

    simulation = archive.data
    print(simulation)
    assert simulation is not None, 'No simulation data in archive'

    assert simulation.model_method, 'No model_method in simulation'

    method = simulation.model_method[0]
    assert hasattr(method, 'numerical_settings'), (
        'ModelMethod missing numerical_settings attribute'
    )
    assert method.numerical_settings is not None, 'numerical_settings is None'
    assert len(method.numerical_settings) > 0, 'numerical_settings is empty'

    # Check for Pseudopotentials
    pseudopotentials = [
        ns for ns in method.numerical_settings if type(ns).__name__ == 'Pseudopotential'
    ]
    assert len(pseudopotentials) > 0, 'No Pseudopotential objects in numerical_settings'

    # Verify OUTCAR-specific fields are populated
    for pp in pseudopotentials:
        assert pp.name is not None
        assert pp.n_valence_electrons is not None
