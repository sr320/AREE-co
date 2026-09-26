# Technical Architecture

The CLI is implemented in Python under `src/aree`.

- `validation`: JSON Schema validation.
- `intake`: registry ingestion.
- `harmonize`: processed-result conversion and identifier mapping.
- `meta_analysis`: random-effects synthesis.
- `prioritize`: candidate scoring.
- `reporting`: evidence cards and demo reports.
- `raw`: shared helpers for raw-data reanalysis scripts: ENA run manifests, checksum-verified resumable FASTQ downloads, sample sheets, RefSeq GFF parsing, and DESeq2-to-evidence export.

Study-specific scripts in `scripts/` keep each study's contrasts, QC gates and report text, and import the shared machinery from `aree.raw`, so they work both as `python scripts/NAME.py` and when imported.

Raw RNA-seq reanalysis runs through the Nextflow pipeline in `workflows/rnaseq/` (see its README). Methylation, proteomics and metabolomics studies currently enter through processed-result harmonization; raw-data pipelines for them are on the roadmap.

