import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger

from tests.parsers.common import SimulationParserPipelineTestSuite


@pytest.mark.large_fixture
class TestCMinimalPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'c_minimal_archive'
    expected_program_name = 'exciting'

    @pytest.fixture(scope='class')
    def normalized_archive(self, archive):
        normalized = EntryArchive.m_from_dict(archive.m_to_dict())
        normalized.metadata = EntryMetadata()
        normalize_all(normalized, logger=get_logger(__name__))
        return normalized


    @pytest.mark.pipeline
    def test_xc_functional_normalization(self, normalized_archive):
        method = normalized_archive.data.model_method[0]
        assert method is not None
        assert method.jacobs_ladder == 'GGA'
        xc = method.xc
        assert xc is not None
        assert xc.functional_key == 'PBESOL'
        assert len(xc.components) == 2
        
        components = {
            component.canonical_label: component for component in xc.components
        }
        assert set(components) == {
            'XC_GGA_C_PBE_SOL',
            'XC_GGA_X_PBE_SOL',
        }

        correlation = components['XC_GGA_C_PBE_SOL']
        assert correlation.libxc_id == 133
        assert correlation.display_name == ('Perdew, Burke & Ernzerhof SOL')
        assert correlation.family == 'GGA'
        assert correlation.kind == 'correlation'

        exchange = components['XC_GGA_X_PBE_SOL']
        assert exchange.libxc_id == 116
        assert exchange.display_name == ('Perdew, Burke & Ernzerhof SOL')
        assert exchange.family == 'GGA'
        assert exchange.kind == 'exchange'


    @pytest.mark.pipeline
    def test_basis_set_container_normalization(self, normalized_archive):
        method = normalized_archive.data.model_method[0]
        assert method is not None
        containers = [
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'BasisSetContainer'
        ]
        assert len(containers) == 1
        assert containers[0].name == 'APW-like'
        components = [
            component
            for component in containers[0].basis_set_components
            if component.m_def.name == 'APWPlaneWaveBasisSet'
        ]
        assert len(components) == 1
        assert components[0].label == "g-mesh"
    
    def test_model_system_normalization(self, normalized_archive):
        system = normalized_archive.data.model_system[0]
        assert system is not None
        assert system.is_representative is True
        assert system.branch_depth == 0
        assert system.composition_formula == 'C(2)'
        assert system.particle_states[0].atomic_number == 6
       

