# Testing Plan for PPCutoff Nested Queries

## Current Status

`PPCutoff.cutoffs` subsection (line 1193 in `numerical_settings.py`) is **NOT configured as nested** - it has no `a_elasticsearch` annotation. This means:

- Elasticsearch will index it as a **flattened array**
- Queries will lose correlation between `cutoff_kind`, `cutoff_role`, and `value` at the same array position
- Cannot query "wavefunction cutoffs with role=recommended AND value > 300 eV" reliably

## Option 1: Test Current Flattened Behavior (No Changes)

### Via API

```bash
# Query any PP with cutoff value > 500 eV (≈0.8e-17 J)
curl -X POST "http://localhost:8000/v1/entries/query" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "results.method.numerical_settings.cutoffs.value": {"gte": 8e-18}
    }
  }'
```

### Via GUI

1. Navigate to search interface
2. Use filter menu to select `Numerical Settings > Cutoffs > Value`
3. Set range filter

### Limitation

Cannot correlate multiple fields (e.g., "wavefunction cutoffs > 500 eV")

## Option 2: Add Nested Support + Test (Recommended)

### Changes Required

1. Add `a_elasticsearch=Elasticsearch(nested=True)` to `cutoffs` SubSection
2. Re-create Elasticsearch index
3. Re-index VASP test data

### Testing Workflow

1. Update schema annotation
2. Parse VASP OUTCAR to populate PPCutoff data
3. Index the entry
4. Query with nested syntax to verify correlation is preserved

### Example Nested Query

```bash
curl -X POST "http://localhost:8000/v1/entries/query" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "results.method.numerical_settings.cutoffs": {
        "all": [
          {"cutoff_kind": "wavefunction"},
          {"cutoff_role": "recommended"},
          {"value": {"gte": 8e-18}}
        ]
      }
    }
  }'
```

## Recommendation

### Start with Option 1

Verify:
- Parser correctly populates PPCutoff subsections
- Data appears in Elasticsearch
- Basic queries work

### Then Implement Option 2

If you need to:
- Query correlated fields (e.g., "wavefunction cutoffs with recommended role")
- Support advanced filtering in GUI
- Enable precise scientific queries

## Testing Commands

### Check Elasticsearch Mapping

```bash
docker exec nomad_elastic curl -X GET \
  "localhost:9200/nomad_entries_v1/_mapping?pretty" | \
  grep -A 10 "cutoffs"
```

### Direct Elasticsearch Query (Flattened)

```bash
docker exec nomad_elastic curl -X GET \
  "localhost:9200/nomad_entries_v1/_search?pretty" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "range": {
        "results.method.numerical_settings.cutoffs.value": {
          "gte": 8e-18
        }
      }
    }
  }'
```

### Direct Elasticsearch Query (Nested)

```bash
docker exec nomad_elastic curl -X GET \
  "localhost:9200/nomad_entries_v1/_search?pretty" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "nested": {
        "path": "results.method.numerical_settings.cutoffs",
        "query": {
          "bool": {
            "must": [
              {"term": {"results.method.numerical_settings.cutoffs.cutoff_kind": "wavefunction"}},
              {"term": {"results.method.numerical_settings.cutoffs.cutoff_role": "recommended"}},
              {"range": {"results.method.numerical_settings.cutoffs.value": {"gte": 8e-18}}}
            ]
          }
        }
      }
    }
  }'
```

## Implementation Steps for Nested Support

### 1. Update Schema

In `packages/nomad-schema-plugins-simulations/src/nomad_simulations/schema_packages/numerical_settings.py`:

```python
from nomad.metainfo import Elasticsearch

# Line 1193 - replace:
cutoffs = SubSection(sub_section=PPCutoff.m_def, repeats=True)

# With:
cutoffs = SubSection(
    sub_section=PPCutoff.m_def,
    repeats=True,
    a_elasticsearch=Elasticsearch(nested=True)
)
```

### 2. Re-create Index

```bash
# Delete existing index
docker exec nomad_elastic curl -X DELETE "localhost:9200/nomad_entries_v1"

# Restart NOMAD to recreate with new mapping
docker compose restart nomad
```

### 3. Re-index Test Data

```bash
uv run nomad parse /home/nathan/Downloads/vasp_example/size_2
```

### 4. Verify Nested Mapping

```bash
docker exec nomad_elastic curl -X GET \
  "localhost:9200/nomad_entries_v1/_mapping?pretty" | \
  grep -A 20 "cutoffs"

# Should show: "type": "nested"
```

## References

- NOMAD Elasticsearch extension: `packages/nomad-FAIR/nomad/metainfo/elasticsearch_extension.py`
- NOMAD search API: `packages/nomad-FAIR/nomad/search.py`
- Your research: `texts/elasticsearch-nested-vs-flattened-arrays.md`
