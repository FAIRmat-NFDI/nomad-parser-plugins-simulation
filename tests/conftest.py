import pytest

# Tests that fail on `test-oasis` because the legacy parsers run first and
# populate the archive before the new parsers.
LEGACY_FIRST_SKIPS = {
    'tests/parsers/fhiaims/test_integration.py': [
        'TestFeScfSpinPolarized::test_method_and_eigenvalue_mapping',
        'TestH2OMolecularDynamics::test_trajectory_mapping',
        'TestSiGeometryOptimization::test_archive_contract',
        'TestSiliconGWBands::test_atoms_k_mesh_band_values_and_occupations',
    ],
    'tests/parsers/h5md/test_integration.py': ['TestH5MDArchive'],
    'tests/parsers/h5md/test_pipeline.py': ['TestH5MDPipeline'],
    'tests/parsers/quantumespresso/test_integration.py': ['TestGIPAWNMRXMLArchive'],
    'tests/parsers/vasp/test_integration.py': [
        'TestAgacRelaxOutcarArchive::test_workflow_convergence',
        'TestAgacRelaxOutcarArchive::test_workflow_has_convergence_targets',
        'TestAgacRelaxVasprunArchive::test_workflow_convergence',
        'TestAgacRelaxVasprunArchive::test_workflow_has_convergence_targets',
    ],
}


def pytest_collection_modifyitems(config, items):
    marker = pytest.mark.skip(reason='Fails with the legacy parsers executed first.')
    for item in items:
        path, _, name = item.nodeid.partition('::')
        name = name.split('[')[0]
        for prefix in LEGACY_FIRST_SKIPS.get(path, []):
            if name == prefix or name.startswith(f'{prefix}::'):
                item.add_marker(marker)
