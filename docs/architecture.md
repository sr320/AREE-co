# Technical Architecture

The CLI is implemented in Python under `src/aree`.

- `validation`: JSON Schema validation.
- `intake`: registry ingestion.
- `harmonize`: processed-result conversion and identifier mapping.
- `meta_analysis`: random-effects synthesis.
- `prioritize`: candidate scoring.
- `reporting`: evidence cards and demo reports.

Raw-data reanalysis is represented by modular Nextflow scaffolds in `workflows/`.

