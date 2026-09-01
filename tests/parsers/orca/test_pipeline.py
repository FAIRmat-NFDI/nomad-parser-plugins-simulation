import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger

from tests.parsers.common import SimulationParserPipelineTestSuite


class TestRIMP2WaterPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'ri_mp2_water_archive'
    expected_program_name = 'ORCA'

    @pytest.mark.pipeline
    def test_frontier_gap_derived(self, archive):
        # ORCA reports no HOMO-LUMO gap directly; it only prints the ORBITAL
        # ENERGIES block. The schema derives the frontier gap in
        # `MolecularOrbitals.normalize()`, reached by the `MetainfoNormalizer` pass
        # (not `Simulation.normalize()`), from `value` and `occupations` of the
        # canonical orbitals. Drive the full normalization the way processing does
        # and assert the gap end-to-end. RI_MP2_water has a clean 2 -> 0 boundary.
        normalized = EntryArchive.m_from_dict(archive.m_to_dict())
        normalized.metadata = EntryMetadata()
        normalize_all(normalized, logger=get_logger(__name__))

        molecular_orbitals = normalized.data.outputs[0].molecular_orbitals[0]
        homo = molecular_orbitals.homo_normalized
        lumo = molecular_orbitals.lumo_normalized
        gap = molecular_orbitals.homo_lumo_gap_normalized
        assert homo is not None and lumo is not None and gap is not None
        assert gap.magnitude >= 0
        assert gap.magnitude == pytest.approx((lumo - homo).magnitude)
