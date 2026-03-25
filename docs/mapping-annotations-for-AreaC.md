# Practical Guide to Mapping Annotations for Parser dev

This doc explains how to populate the NOMAD simulations datamodel using mapping annotations. It shows what they are, how to write them, and how to handle nested and repeating sections. 

---

## What a mapping annotation is

- A small note on a schema field that says: *for source `<KEY>`, read value at `<PATH>` and put it here.*
- Stored on any section/quantity: `m_annotations['mapping'][<KEY>] = MapperAnnotation(mapper=...)`.
- `<KEY>` names the source (e.g., `OUT_KEY`, `HDF5_KEY`, `INFO_KEY`, `BAND_KEY`).
- Paths use JMESPath style. Relative paths start with `.` (from the parent’s mapped node). `@` anchors at the source root.
- The machinery (`MetainfoParser`) walks the schema, follows these notes, instantiates sections/subsections, and fills values.

## What happens under the hood

1) A source parser (`XMLParser`, `HDF5Parser`, `TextParser`, etc.) reads a file into a Python dict.
2) Set `annotation_key` on a `MetainfoParser` to the `<KEY>` that matches the annotations.
3) Call `convert(source_parser, target_parser)`. One key → one `convert`. Multiple keys/sources → run `convert` once per key after switching `annotation_key` each time.
4) The target archive is populated; repeating subsections are created automatically from lists.

## Advantages of using mapping

- Declarative and readable: data movement is described alongside the schema, not scattered in code.
- Reusable across sources: same schema, different annotation keys; swap sources by switching `annotation_key`.
- Automatic instantiation: nested and repeating subsections are created based on list structure; no manual wiring.
- Composable transforms: reshape, filter, and unit application via small helper functions in the parser.
- Easier maintenance: when input formats change, update paths in one place (annotations) instead of parser code.
- Parallel branches: ingest multiple files/sources into the same archive by running `convert` per key.

## Path cheat sheet

- Relative from parent: `.field`, `.field[0]`, `.list[*].value`
- Absolute from root: `field.subfield` (no leading dot)
- Root anchor: `@`
- Filter: `.items[?"@name"=="version"] | [0].__value`
- Attributes (XML/HDF5): keys often prefixed with `@`; values with `__value`

## When to annotate `m_def` vs. the attribute

- Attribute-level (on the quantity or subsection attribute):
  - Use for a specific, local mapping. This is the default for quantities. Example: `Program.version` with `.program_version`.
  - It overrides `m_def` if both exist.
- `m_def` (on the section definition):
  - Use for a default that should work anywhere that section appears. Example: `Simulation.m_def` with `@` to anchor the whole tree.
  - Use to give a subsection its own base path so children can resolve relative paths. Example: `DFT.m_def` with `.dft`, so all quantities inside DFT can start from `.dft` even if the parent doesn’t have a special annotation.
  - Use for abstract/base sections when you need a mapping that applies to all their concrete inheritors.
- Lookup order: attribute annotation → `m_def` annotation → inheriting section definitions (for abstract subsections). If you want a local override, put it on the attribute; if you want a reusable default or a base path, put it on `m_def`.

## Core patterns 

### 1) Root anchor
Anchor the target root section to the **source tree root** (`@` = top of source dict/XML/HDF5).
```python
add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
```

### 2) Nested quantities
Relative to the parent section’s mapped node.
```python
add_mapping_annotation(general.Program.name, OUT_KEY, '.program')
add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')
```

### 3) Single nested subsection
Annotate the subsection’s section definition.
```python
add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.dft')
```

### 4) Repeating subsections (lists)
Point to a list; one subsection is created per element. For repeating quantities or array-like quantities, map to a list/array; it is set directly.
```python
add_mapping_annotation(model_system.AtomsState.m_def, OUT_KEY, '.atoms')
add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.label')
add_mapping_annotation(model_system.AtomsState.atomic_number, OUT_KEY, '.number')
add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')  # shape [n_atoms, 3]
```

### 5) Multiple sources into the same schema
Use different keys; run `convert` once per key.
```python
# Octopus uses OUT_KEY + INFO_KEY + EIGENVALUES_KEY on the same schema
# Exciting uses INFO_KEY, INPUT_XML_KEY, EIGVAL_KEY, BANDSTRUCTURE_XML_KEY, DOS_XML_KEY
# Wannier90 uses WOUT_KEY, WIN_KEY, BAND_KEY, WHR_KEY, DOS_KEY
```
Set `annotation_key = <KEY>` before each `convert`.

### 6) Transformations (reshape, filter, apply units)
Use a tuple `(function_name, [paths], {kwargs})` and implement the function on the parser.
```python
add_mapping_annotation(
    outputs.TotalForce.value,
    SECONDARY_KEY,
    ('reshape_forces', ['.raw_forces'], {'unit': 'eV / angstrom'})
)
```


### 7) Units or subfield search
You can pass `unit='eV'` or `search='.path.to.value'` in the annotation (see Wannier90 eigenvalue/DOS paths).

## Starter template

File: `packages/nomad_simulation_parsers/schema_packages/template.py`
```python
from nomad_simulation_parsers.schema_packages import general, model_system, outputs
from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

PRIMARY_KEY = 'primary'      # main source (text/JSON/XML)
SECONDARY_KEY = 'secondary'  # auxiliary source

# Section anchor
add_mapping_annotation(general.Simulation.m_def, PRIMARY_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, SECONDARY_KEY, '@')

# Nested quantities
add_mapping_annotation(general.Program.name, PRIMARY_KEY, '.program.name')
add_mapping_annotation(general.Program.version, PRIMARY_KEY, '.program.version')

# Single nested subsection
add_mapping_annotation(model_system.ModelSystem.m_def, PRIMARY_KEY, '.structure')
add_mapping_annotation(model_system.ModelSystem.lattice_vectors, PRIMARY_KEY, '.lattice')
add_mapping_annotation(model_system.ModelSystem.periodic_boundary_conditions, PRIMARY_KEY, '.pbc')

# Repeating subsections
add_mapping_annotation(model_system.AtomsState.m_def, PRIMARY_KEY, '.atoms')
add_mapping_annotation(model_system.AtomsState.chemical_symbol, PRIMARY_KEY, '.symbol')
add_mapping_annotation(model_system.AtomsState.atomic_number, PRIMARY_KEY, '.Z')
add_mapping_annotation(model_system.ModelSystem.positions, PRIMARY_KEY, '.positions')

# Repeating outputs
add_mapping_annotation(outputs.Outputs.total_energies, PRIMARY_KEY, '.energies')
add_mapping_annotation(outputs.TotalEnergy.value, PRIMARY_KEY, '.value')
add_mapping_annotation(outputs.TotalEnergy.name, PRIMARY_KEY, '.label')

# Transformation example on secondary source
add_mapping_annotation(
    outputs.TotalForce.value,
    SECONDARY_KEY,
    ('reshape_forces', ['.raw_forces'], {'unit': 'eV / angstrom'}),
)
```

File: `packages/nomad_simulation_parsers/parsers/template/parser.py`
```python
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad_simulations.schema_packages.general import Simulation, Program
from nomad_simulation_parsers.schema_packages import template  # ensure annotations load

class TemplateTextParser(TextParser):
    pass  # plug in your actual text parser; must expose ._results

class TemplateParser:
    def __init__(self):
        self._source_primary = TemplateTextParser()
        self._source_secondary = TemplateTextParser()  # swap for XML/HDF5 parsers as needed
        self._target = MetainfoParser()
        self.archive = None

    @staticmethod
    def reshape_forces(raw, unit=None):
        import numpy as np
        if raw is None:
            return None
        arr = np.array(raw).reshape(-1, 3)
        return arr  # unit handling can be added if needed

    def parse(self, mainfile):
        self._source_primary.filepath = mainfile
        self._source_secondary.filepath = mainfile.replace('.out', '.aux')  # adjust

        self.archive = Simulation(program=Program(name='TEMPLATE'))
        self._target.data_object = self.archive

        self._target.annotation_key = template.PRIMARY_KEY
        self._source_primary.convert(self._target)

        self._target.annotation_key = template.SECONDARY_KEY
        self._source_secondary.convert(self._target)

        return self.archive
```

## How to add a new mapping (step-by-step)

1) Pick the target field in `nomad-simulations` (quantity or subsection).
2) Choose a source key (`OUT_KEY`, `HDF5_KEY`, etc.).
3) Write the JMESPath that matches your source dict. Use `.` to start from the parent node.
4) Add the annotation with `add_mapping_annotation(target, KEY, path_or_tuple)`.
5) Repeats: if the subsection is `repeats=True`, just point to a list path; one instance per list element.
6) If combining multiple sources, run `convert` for each key (set `annotation_key` before calling).
7) Test on a small file, inspect the archive (e.g., print `archive.data`), and adjust paths.

## Pitfalls to avoid

- Forgetting the leading `.` on relative paths (they become absolute and miss).
- Abstract subsections: put annotations on concrete inheriting sections or rely on inheriting annotations.
- Reusing a schema with different sources: clear old mappings first with `remove_mapping_annotations`.
- Shape issues: for reshaping or unit conversion, use a transformer tuple and implement the helper.
