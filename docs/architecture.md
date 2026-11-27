# Technical Architecture

The CLI is implemented in Python under `src/aree`.

- `validation`: JSON Schema validation.
- `intake`: registry ingestion.
- `harmonize`: processed-result conversion and identifier mapping.
- `meta_analysis`: random-effects synthesis.
- `prioritize`: candidate scoring.
- `reporting`: evidence cards, pure-Python SVG forest plots (`forest_plot.py`), and demo reports.
- `export`: relational SQLite, ZSTD Parquet + Frictionless `datapackage.json`, and Schema.org JSON-LD catalog exporters (`database.py`).
- `raw`: shared helpers for raw-data reanalysis scripts: ENA run manifests, checksum-verified resumable FASTQ downloads, NCBI SRA verification and streaming conversion, versioned multi-species RefSeq bundle fetching (`reference.py`), multi-assay evidence validation/export (`assay_evidence.py`), sample sheets, RefSeq GFF parsing, DESeq2-to-evidence export, and annotation table merging.

Study-specific scripts in `scripts/` keep each study's contrasts, QC gates and report text, and import the shared machinery from `aree.raw`, so they work both as `python scripts/NAME.py` and when imported.

Reusable Nextflow pipelines live in `workflows/rnaseq/` (paired/single-end, blocking covariates, replicate collapsing, 3' tag-seq, in-pipeline `--reference_accession` resolution via `modules/fetch_reference.nf`), `workflows/methylation/` (Bismark / CpG M-value `limma`), `workflows/proteomics/` (LFQ/TMT `limma`), and `workflows/metabolomics/` (PQN + `limma`/Wilcoxon), each using `modules/validate_evidence.nf` and container definitions in `containers/`.

