# Technical Architecture

The CLI is implemented in Python under `src/aree`.

- `validation`: JSON Schema validation.
- `intake`: registry ingestion.
- `harmonize`: processed-result conversion and identifier mapping.
- `meta_analysis`: random-effects synthesis.
- `prioritize`: candidate scoring.
- `reporting`: evidence cards and demo reports.
- `raw`: shared helpers for raw-data reanalysis scripts: ENA run manifests, checksum-verified resumable FASTQ downloads, NCBI SRA verification and streaming conversion, sample sheets, RefSeq GFF parsing, DESeq2-to-evidence export, and annotation table merging.

Study-specific scripts in `scripts/` keep each study's contrasts, QC gates and report text, and import the shared machinery from `aree.raw`, so they work both as `python scripts/NAME.py` and when imported.

The reusable Nextflow pipeline in `workflows/rnaseq/` supports paired-end and single-end RNA-seq reanalysis, blocking covariates, technical-replicate collapsing, and 3' tag-seq (`--tag_seq true`) (see its README). Committed public-study analyses also preserve per-study shell wrappers in `scripts/` for provenance. Methylation, proteomics and metabolomics studies currently enter through processed-result harmonization; raw-data pipelines for them are on the roadmap.

