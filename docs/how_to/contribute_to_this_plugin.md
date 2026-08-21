# Contribute to This Plugin

Before opening a pull request, run the relevant tests and include any generated
files that are affected by your changes.

## Update the parser mapping report

The parser mapping report is a committed snapshot of the file-parser quantities
and their archive mappings. Its generation is part of CI. Whenever a pull
request contains source-code changes to be merged into `develop`, regenerate
the report and commit the resulting update as part of that pull request:

```sh
nomad-sim-parser mapping-report \
  --override docs/reference/parser_mapping_report_overrides.yaml
```

This updates
`docs/reference/parser_mapping_report.md`. Review the generated diff, including
any intentionally unmapped quantities, before submitting the pull request.
