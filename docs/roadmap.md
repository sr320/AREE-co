# Roadmap

## MVP

- Study schema and validation.
- Controlled vocabularies.
- Synthetic *C. gigas* demo studies.
- Processed-result harmonization.
- Random-effects meta-analysis.
- Candidate scoring.
- Evidence cards.
- Streamlit browser and Markdown reports.

## Available real-study workflows

- Verified public Pacific oyster accessions and documented study designs.
- Raw RNA-seq reanalysis and harmonized evidence, with a synthetic end-to-end Nextflow check in CI.
- Study-specific QuantSeq 3' Tag-seq reanalysis with counts exported without a transcript-length offset.
- Tank-aware mixed-model RNA-seq reanalysis (`variancePartition::dream` with Satterthwaite degrees of freedom) for `PRJNA735889`.
- Real-study findings, comparability-reviewed cross-study synthesis (`config/synthesis_strata.yaml`), evidence cards, provenance and downloads on the primary website; synthetic results in a secondary demo section.
- Covariance-aware DL and REML meta-analysis (`--shared-control-rho`) cross-checked against R `metafor`, with candidate scoring weight sensitivity audits.
- Multi-species RefSeq reference catalog (`config/references.yaml`, `aree fetch-reference`, `modules/fetch_reference.nf`) and OrthoFinder 3.1.5 8-bivalve ortholog mapping releases (`data/mappings/orthologs/`).
- Containerized multi-assay Nextflow pipelines (`workflows/methylation/`, `workflows/proteomics/`, `workflows/metabolomics/`, `containers/`) with post-export evidence validation (`modules/validate_evidence.nf`).
- Pure-Python SVG forest plots embedded in evidence cards, Quarto cross-study synthesis pages, and the Streamlit app (`src/aree/reporting/forest_plot.py`).
- Database export targets (`aree export-db --format sqlite|parquet|jsonld`) with Frictionless Data Package and Schema.org JSON-LD support (`reports/exports/`).

## Production Next Steps

- Expand the verified public-study collection across bivalve species, prioritizing measured resilience outcomes, multi-omics integration, and independent replication.
