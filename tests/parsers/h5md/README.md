# H5MD parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| H5MD mainfile recognition | not applicable | not applicable | representative HDF5 fixture |
| H5MD header and program metadata | not applicable | `test_mapping.py` | `test_h5md_archive_contract` |
| Particle trajectory and units | not applicable | `test_mapping.py` | `test_h5md_archive_contract` |
| Topology, particle identity, and subsystem hierarchy | not applicable | `test_mapping.py` | `test_h5md_archive_contract` |
| Configurational outputs and contributions | not applicable | `test_mapping.py` | `test_h5md_archive_contract` |
| Workflow method and MD results | not applicable | `test_mapping.py` | `test_h5md_workflow_contract` |
| Ensemble and correlation-function results | not applicable | `test_mapping.py` | `test_h5md_workflow_contract` |
| Invalid or inconsistent trajectory data | not applicable | `test_mapping.py` | not applicable |
| NOMAD normalization and identity placement | not applicable | not applicable | `test_pipeline.py` |

H5MD uses the standard HDF5 reader from `mapping_parser`; its parser-specific
unit tests therefore target only the source-to-metainfo mapping callbacks. The
complete fixture remains an integration test because it exercises the topology,
trajectory, observables, and workflow mappings together.
