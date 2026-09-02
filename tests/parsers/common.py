"""Reusable contracts for public simulation parser archives.

The assertions in this module deliberately stop at parser-independent NOMAD
invariants. Source-format recognition and scientific values belong in each
parser's own tests.
"""

import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger


class _SimulationParserSuite:
    """Configuration and fixtures shared by parser contract suites."""

    archive_fixture: str
    expected_program_name: str

    @pytest.fixture(scope='class')
    def archive(self, request) -> EntryArchive:
        return request.getfixturevalue(self.archive_fixture)


class SimulationParserTestSuite(_SimulationParserSuite):
    """Common public-archive integration contracts for a parser.

    Subclasses only need to name a module-scoped archive fixture and the
    expected program name. The class name intentionally does not start with
    ``Test`` so pytest only collects configured subclasses.
    """

    required_simulation_sections = ('model_method', 'model_system', 'outputs')
    require_lattice_vectors = False
    require_periodic_boundary_conditions = False

    @pytest.mark.integration
    def test_archive_is_valid(self, archive):
        errors, warnings = archive.m_validate()

        assert errors == []
        assert warnings == []

    @pytest.mark.integration
    def test_archive_has_required_sections(self, archive):
        simulation = archive.data

        assert simulation is not None
        assert simulation.program is not None
        assert simulation.program.name == self.expected_program_name
        assert simulation.model_system
        for section_name in self.required_simulation_sections:
            assert getattr(simulation, section_name), (
                f'missing required Simulation.{section_name} section'
            )

    @pytest.mark.integration
    def test_archive_serialization_round_trip(self, archive):
        restored = EntryArchive.m_from_dict(archive.m_to_dict())
        errors, warnings = restored.m_validate()

        assert errors == []
        assert warnings == []
        assert restored.data.program.name == self.expected_program_name
        assert len(restored.data.model_method) == len(archive.data.model_method)
        assert len(restored.data.outputs) == len(archive.data.outputs)

    @pytest.mark.integration
    def test_archive_has_model_systems(self, archive):
        assert archive.data is not None
        assert archive.data.model_system

    @pytest.mark.integration
    def test_representative_system_is_complete(self, archive):
        simulation = archive.data
        representative = next(
            (
                system
                for system in simulation.model_system
                if getattr(system, 'is_representative', False)
            ),
            simulation.model_system[0],
        )

        assert representative.positions is not None
        assert representative.particle_states
        assert all(
            state.chemical_symbol is not None
            for state in representative.particle_states
        )
        if self.require_periodic_boundary_conditions:
            assert representative.periodic_boundary_conditions is not None
        if self.require_lattice_vectors:
            assert representative.lattice_vectors is not None
        if representative.periodic_boundary_conditions is not None and any(
            representative.periodic_boundary_conditions
        ):
            assert representative.lattice_vectors is not None

    @pytest.mark.integration
    def test_model_system_serialization_round_trip(self, archive):
        restored = EntryArchive.m_from_dict(archive.m_to_dict())

        assert len(restored.data.model_system) == len(archive.data.model_system)


class WorkflowTestSuite(_SimulationParserSuite):
    """Common workflow contracts for public parser archives."""

    workflow_name: str = 'SimulationWorkflow'

    @pytest.mark.integration
    def test_archive_has_workflow(self, archive):
        assert archive.workflow2 is not None
        assert archive.workflow2.m_def.name == self.workflow_name

    @pytest.mark.integration
    def test_workflow_serialization_round_trip(self, archive):
        restored = EntryArchive.m_from_dict(archive.m_to_dict())

        assert restored.workflow2.m_def.name == archive.workflow2.m_def.name


class SimulationParserPipelineTestSuite(_SimulationParserSuite):
    """Common NOMAD normalization compatibility contract for a parser."""

    @pytest.mark.pipeline
    def test_normalizer_consumes_archive(self, archive):
        normalized = EntryArchive.m_from_dict(archive.m_to_dict())
        normalized.metadata = EntryMetadata()

        normalize_all(normalized, logger=get_logger(__name__))

        representative_index = normalized.data.representative_system_index
        assert representative_index is not None
        assert 0 <= representative_index < len(normalized.data.model_system)
        assert normalized.results is not None
        assert normalized.results.properties is not None
