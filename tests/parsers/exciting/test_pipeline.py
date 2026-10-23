import pytest

from tests.parsers.common import SimulationParserPipelineTestSuite


@pytest.mark.large_fixture
class TestCMinimalPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'c_minimal_archive'
    expected_program_name = 'exciting'


    @pytest.mark.pipeline
    def test_xc_functional_normalization(self, archive):
        from nomad.client import normalize_all
        from nomad.datamodel import EntryArchive, EntryMetadata
        from nomad.utils import get_logger

        normalized = EntryArchive.m_from_dict(archive.m_to_dict())
        normalized.metadata = EntryMetadata()

        normalize_all(normalized, logger=get_logger(__name__))

        method = normalized.data.model_method[0]
        xc = method.xc

        assert method.jacobs_ladder == 'GGA'
        assert xc is not None
        assert xc.functional_key == 'PBESOL'
        assert len(xc.components) == 2

        components = {
        component.canonical_label: component
        for component in xc.components
        }

        assert set(components) == {
        'XC_GGA_C_PBE_SOL',
        'XC_GGA_X_PBE_SOL',
        }

        correlation = components['XC_GGA_C_PBE_SOL']
        assert correlation.libxc_id == 133
        assert correlation.display_name == (
        'Perdew, Burke & Ernzerhof SOL'
        )
        assert correlation.family == 'GGA'
        assert correlation.kind == 'correlation'

        exchange = components['XC_GGA_X_PBE_SOL']
        assert exchange.libxc_id == 116
        assert exchange.display_name == (
        'Perdew, Burke & Ernzerhof SOL'
        )
        assert exchange.family == 'GGA'
        assert exchange.kind == 'exchange'

