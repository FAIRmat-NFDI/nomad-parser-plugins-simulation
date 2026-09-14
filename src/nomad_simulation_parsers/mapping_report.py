"""Generate file-parser to archive-mapper dependency reports."""

import argparse
import importlib
from dataclasses import dataclass
from pathlib import Path

import yaml
from nomad_file_parser.mapping_parser import MetainfoParser
from nomad_simulations.schema_packages.general import Simulation

from nomad_simulation_parsers.parsers.utils.general import create_mapping_table


@dataclass(frozen=True)
class ParserReportSpec:
    """A file-parser class and the archive mapping key it feeds."""

    name: str
    file_parser_module: str
    file_parser_class: str
    schema_module: str
    annotation_key: str


# Keep the associations explicit: a parser can expose several files, each with
# its own mapping annotation key.  Additional parsers can be added here without
# changing the report generation logic.
PARSER_REPORT_SPECS = (
    ParserReportSpec(
        'ABINIT / OUT',
        'nomad_simulation_parsers.parsers.abinit.file_parser',
        'AbinitOutParser',
        'nomad_simulation_parsers.schema_packages.abinit',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'AMS / OUT',
        'nomad_simulation_parsers.parsers.ams.file_parser',
        'OutParser',
        'nomad_simulation_parsers.schema_packages.ams',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Crystal / OUT',
        'nomad_simulation_parsers.parsers.crystal.file_parser',
        'OutputParser',
        'nomad_simulation_parsers.schema_packages.crystal',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Crystal / F25',
        'nomad_simulation_parsers.parsers.crystal.file_parser',
        'F25Parser',
        'nomad_simulation_parsers.schema_packages.crystal',
        'F25_KEY',
    ),
    ParserReportSpec(
        'Exciting / INFO.OUT',
        'nomad_simulation_parsers.parsers.exciting.info_parser',
        'InfoFileParser',
        'nomad_simulation_parsers.schema_packages.exciting',
        'INFO_KEY',
    ),
    ParserReportSpec(
        'Exciting / EIGVAL.OUT',
        'nomad_simulation_parsers.parsers.exciting.eigval_parser',
        'EigvalFileParser',
        'nomad_simulation_parsers.schema_packages.exciting',
        'EIGVAL_KEY',
    ),
    ParserReportSpec(
        'FHI-aims / output',
        'nomad_simulation_parsers.parsers.fhiaims.out_parser',
        'FHIAimsOutFileParser',
        'nomad_simulation_parsers.schema_packages.fhiaims',
        'TEXT_KEY',
    ),
    ParserReportSpec(
        'GPAW / GPW',
        'nomad_simulation_parsers.parsers.gpaw.gpw_parser',
        'GPWFileParser',
        'nomad_simulation_parsers.schema_packages.gpaw',
        'GPW_KEY',
    ),
    ParserReportSpec(
        'LOBSTER / OUT',
        'nomad_simulation_parsers.parsers.lobster.file_parser',
        'OutParser',
        'nomad_simulation_parsers.schema_packages.lobster',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'LOBSTER / COXPCAR',
        'nomad_simulation_parsers.parsers.lobster.file_parser',
        'COXPCARParser',
        'nomad_simulation_parsers.schema_packages.lobster',
        'COXPCAR_KEY',
    ),
    ParserReportSpec(
        'LOBSTER / CHARGE',
        'nomad_simulation_parsers.parsers.lobster.file_parser',
        'CHARGEParser',
        'nomad_simulation_parsers.schema_packages.lobster',
        'CHARGE_KEY',
    ),
    ParserReportSpec(
        'LOBSTER / ICOHPLIST',
        'nomad_simulation_parsers.parsers.lobster.file_parser',
        'ICOXPLISTParser',
        'nomad_simulation_parsers.schema_packages.lobster',
        'ICOXPLIST_KEY',
    ),
    ParserReportSpec(
        'Octopus / OUT',
        'nomad_simulation_parsers.parsers.octopus.file_parser',
        'OutParser',
        'nomad_simulation_parsers.schema_packages.octopus',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Octopus / eigenvalues',
        'nomad_simulation_parsers.parsers.octopus.file_parser',
        'EigenvalueParser',
        'nomad_simulation_parsers.schema_packages.octopus',
        'EIGENVALUES_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / OUT',
        'nomad_simulation_parsers.parsers.quantumespresso.file_parser',
        'QuantumEspressoFileParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.common',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / DOS',
        'nomad_simulation_parsers.parsers.quantumespresso.pwscf.file_parser',
        'PWSCFFileParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.common',
        'DOS_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / DOS output',
        'nomad_simulation_parsers.parsers.quantumespresso.pwscf.file_parser',
        'PWSCFDOSTextParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.common',
        'DOS_OUT_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / EPW',
        'nomad_simulation_parsers.parsers.quantumespresso.epw.file_parser',
        'EPWFileParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.common',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / phonon',
        'nomad_simulation_parsers.parsers.quantumespresso.phonon.file_parser',
        'PhononFileParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.phonon',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / XSpectra',
        'nomad_simulation_parsers.parsers.quantumespresso.xspectra.file_parser',
        'XSpectraFileParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.common',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Quantum ESPRESSO / GIPAW',
        'nomad_simulation_parsers.parsers.quantumespresso.gipaw.file_parser',
        'GIPAWFileParser',
        'nomad_simulation_parsers.schema_packages.quantumespresso.gipaw',
        'GIPAW_OUT_KEY',
    ),
    ParserReportSpec(
        'VASP / OUTCAR',
        'nomad_simulation_parsers.parsers.vasp.outcar_parser',
        'OutcarTextParser',
        'nomad_simulation_parsers.schema_packages.vasp',
        'OUTCAR_KEY',
    ),
    ParserReportSpec(
        'VASP / vasprun.xml',
        'nomad_simulation_parsers.parsers.vasp.xml_parser',
        'VasprunParser',
        'nomad_simulation_parsers.schema_packages.vasp',
        'XML_KEY',
    ),
    ParserReportSpec(
        'Yambo / OUT',
        'nomad_simulation_parsers.parsers.yambo.file_parsers',
        'MainfileParser',
        'nomad_simulation_parsers.schema_packages.yambo',
        'OUT_KEY',
    ),
    ParserReportSpec(
        'Yambo / NetCDF',
        'nomad_simulation_parsers.parsers.yambo.file_parsers',
        'NetCDFParser',
        'nomad_simulation_parsers.schema_packages.yambo',
        'NETCDF_KEY',
    ),
    ParserReportSpec(
        'Wannier90 / wout',
        'nomad_simulation_parsers.parsers.wannier90.file_parsers',
        'WOutParser',
        'nomad_simulation_parsers.schema_packages.wannier90',
        'WOUT_KEY',
    ),
    ParserReportSpec(
        'Wannier90 / win',
        'nomad_simulation_parsers.parsers.wannier90.file_parsers',
        'WInParser',
        'nomad_simulation_parsers.schema_packages.wannier90',
        'WIN_KEY',
    ),
    ParserReportSpec(
        'Wannier90 / hr',
        'nomad_simulation_parsers.parsers.wannier90.file_parsers',
        'HrParser',
        'nomad_simulation_parsers.schema_packages.wannier90',
        'WHR_KEY',
    ),
)

# These parsers write to the archive through custom code rather than through
# declarative file-parser quantities and archive-mapper annotations.  Keep them
# visible in the report, but do not manufacture a misleading coverage value.
SKIPPED_PARSERS = (
    'LAMMPS',
    'Phonopy',
)


def _read_overrides(path: str | Path | None) -> dict:
    if path is None:
        return {}
    override_path = Path(path).expanduser()
    if not override_path.exists():
        raise FileNotFoundError(f'Override YAML file not found: {override_path}')
    data = yaml.safe_load(override_path.read_text(encoding='utf-8'))
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError('Override YAML must contain a top-level mapping.')
    overrides = data.get('parsers', data) or {}
    if not isinstance(overrides, dict):
        raise ValueError('Override YAML "parsers" value must be a mapping.')
    return overrides


def _apply_quantity_override(
    row: dict,
    parser_name: str,
    override: object,
) -> None:
    if not isinstance(override, dict):
        raise ValueError(
            f'Override for {parser_name!r}/{row["quantity"]!r} must be a mapping.'
        )
    if 'mapped' in override:
        row['mapped'] = bool(override['mapped'])
    if 'mapping' not in override:
        return

    mapping = override['mapping']
    if isinstance(mapping, str):
        mapping = [mapping]
    if not isinstance(mapping, list) or not all(
        isinstance(value, str) for value in mapping
    ):
        raise ValueError(
            f'Mapping override for {parser_name!r}/{row["quantity"]!r} '
            'must be a string or list of strings.'
        )
    row['mapping'] = mapping


def _apply_overrides(
    rows: dict[str, list[dict]],
    skipped: list[str],
    overrides: dict,
) -> None:
    for parser_name, parser_override in overrides.items():
        if not isinstance(parser_override, dict):
            raise ValueError(f'Override for {parser_name!r} must be a mapping.')
        if parser_override.get('skip'):
            rows.pop(parser_name, None)
            if parser_name not in skipped:
                skipped.append(parser_name)
            continue

        quantity_overrides = parser_override.get('quantities', {})
        if not isinstance(quantity_overrides, dict):
            raise ValueError(
                f'Quantity overrides for {parser_name!r} must be a mapping.'
            )
        for row in rows.get(parser_name, []):
            override = quantity_overrides.get(row['quantity'])
            if override is not None:
                _apply_quantity_override(row, parser_name, override)


def generate_rows(
    override_path: str | Path | None = None,
) -> tuple[dict[str, list[dict]], list[str]]:
    rows = {}
    skipped = list(SKIPPED_PARSERS)
    for spec in PARSER_REPORT_SPECS:
        try:
            file_parser_module = importlib.import_module(spec.file_parser_module)
            file_parser = getattr(file_parser_module, spec.file_parser_class)()
            schema = importlib.import_module(spec.schema_module)
            archive_parser = MetainfoParser(data_object=Simulation())
            archive_parser.annotation_key = getattr(schema, spec.annotation_key)
            function_objects = [file_parser_module]
            parser_module_name = f'{spec.file_parser_module.rsplit(".", 1)[0]}.parser'
            try:
                parser_module = importlib.import_module(parser_module_name)
                function_objects.append(parser_module)
                function_objects.extend(
                    value
                    for value in vars(parser_module).values()
                    if isinstance(value, type)
                )
            except ImportError:
                pass
            rows[spec.name] = create_mapping_table(
                file_parser, archive_parser, function_objects=function_objects
            )
        except Exception:
            skipped.append(spec.name)
    _apply_overrides(rows, skipped, _read_overrides(override_path))
    return rows, skipped


def render_table(rows: list[dict]) -> str:
    if not rows:
        return (
            '**Coverage:** Not available. A reliable coverage report cannot be '
            'generated because this parser does not expose reportable file-parser '
            'quantities.'
        )
    mapped = sum(row['mapped'] for row in rows)
    coverage = mapped / len(rows) * 100 if rows else 0
    lines = [
        f'**Summary:** {mapped} mapped, {len(rows) - mapped} unmapped quantities '
        f'({coverage:.2f}% coverage).',
        '',
        '| File-parser quantity | Status | Archive mapper source |',
        '| --- | --- | --- |',
    ]
    for row in rows:
        source = (
            '<br>'.join(f'`{value}`' for value in row['mapping'])
            if row['mapping']
            else '—'
        )
        status = 'Mapped' if row['mapped'] else 'Unmapped'
        lines.append(f'| `{row["quantity"]}` | {status} | {source} |')
    return '\n'.join(lines)


def render_report(rows: dict[str, list[dict]], skipped: list[str] | None = None) -> str:
    sections = [
        '# Simulation parser mapping report',
        '',
        'This report is generated with `create_mapping_table` by comparing the '
        'quantities declared in the simulation file parsers with the source paths '
        'used by their archive mappers.',
        '',
        'The report is a snapshot of the mapper configuration. To inspect another '
        'file parser, call `create_mapping_table(file_parser, archive_parser)` '
        'from `nomad_simulation_parsers.parsers.utils.general` and serialize the '
        'returned rows as Markdown or another table format.',
        '',
        'Optional manual corrections can be supplied with '
        '`nomad-sim-parser mapping-report --override <file.yaml>`; see '
        '`docs/reference/parser_mapping_report_overrides.yaml` for the format.',
        '',
    ]
    for filename, file_rows in rows.items():
        sections.extend([f'## {filename}', '', render_table(file_rows), ''])
    if skipped:
        sections.extend(
            [
                '## Skipped parser specifications',
                '',
                'The following parsers could not be represented by a reliable '
                'file-parser quantity mapping report because they use custom '
                'archive-writing logic or do not expose reportable quantities. '
                'No coverage value was generated for them:',
                '',
                *[f'- `{item}`' for item in skipped],
                '',
            ]
        )
    return '\n'.join(sections)


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the mapping-report command-line parser."""
    parser = argparse.ArgumentParser(
        prog='nomad-sim-parser mapping-report',
        description='Generate the simulation file-parser mapping report.',
    )
    parser.add_argument(
        '--output',
        default='docs/reference/parser_mapping_report.md',
        help='Output Markdown path (default: %(default)s).',
    )
    parser.add_argument(
        '--override',
        default=None,
        help='Optional YAML file with report overrides.',
    )
    return parser


def run(arguments: list[str] | None = None) -> int:
    """Parse mapping-report arguments and write the report."""
    parser = build_argument_parser()
    args = parser.parse_args(arguments)
    return main(args.output, args.override)


def main(
    output: str | Path | None = None,
    override: str | Path | None = None,
) -> int:
    """Write the mapping report and return a process status code."""
    if output is None:
        return run()
    output_path = Path(output).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows, skipped = generate_rows(override)
    output_path.write_text(render_report(rows, skipped), encoding='utf-8')
    print(f'Updated: {output_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(run())
