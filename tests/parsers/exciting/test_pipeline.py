import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger

from tests.parsers.common import SimulationParserPipelineTestSuite


@pytest.mark.large_fixture
class TestCMinimalLargeDataPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'c_minimal_largedata_archive'
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
    def test_model_system_normalization(self, normalized_archive):
        system = normalized_archive.data.model_system[0]
        assert system is not Nonegi
        assert system.is_representative is True
        assert system.branch_depth == 0
        assert system.composition_formula == 'C(2)'
        assert system.particle_states[0].atomic_number == 6
       
    
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
    def test_basis_set_container_normalization(self, normalized_archive):
        method = normalized_archive.data.model_method[0]
        assert method is not None
        containers = [
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'BasisSetContainer'
        ]
        assert len(containers) == 1
        assert containers[0].name == 'LAPW+lo'
        components = [
            component
            for component in containers[0].basis_set_components
            if component.m_def.name == 'APWPlaneWaveBasisSet'
        ]
        assert len(components) == 1
        assert components[0].label == "g-mesh"

    @pytest.mark.pipeline
    def test_orbital_counts_normalization(self, normalized_archive):
        method = normalized_archive.data.model_method[0]
        assert method is not None

        containers = [
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'BasisSetContainer'
        ]
        assert len(containers) == 1

        components = [
            component
            for component in containers[0].basis_set_components
            if component.m_def.name == 'MuffinTinRegion'
        ]
        assert len(components) == 1

        channel = components[0].l_channels[0]
        assert channel.name == 0
        assert channel.n_orbitals == 2
        assert len(channel.orbitals) == 2

        apw = channel.orbitals[0]
        assert apw.m_def.name == 'APWOrbital'
        assert apw.n_terms == 1
        assert len(apw.energy_parameter) == apw.n_terms
        assert len(apw.differential_order) == apw.n_terms

        local_orbital = channel.orbitals[1]
        assert local_orbital.m_def.name == 'APWLocalOrbital'
        assert local_orbital.n_terms == 2
        assert len(local_orbital.energy_parameter) == local_orbital.n_terms
        assert len(local_orbital.differential_order) == local_orbital.n_terms

