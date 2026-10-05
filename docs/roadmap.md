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
- Real-study findings, evidence cards, provenance and downloads on the primary website; synthetic results in a secondary demo section.
- Shared-control replication flags and explicit limitations for the exploratory tank-dependent analysis.

## Production Next Steps

- Expand the verified public-study collection, prioritizing measured resilience outcomes and independent replication.
- Reanalyse PRJNA735889 with a tank-aware model before using its nominal uncertainty for inference.
- Review contrast comparability and model shared-sample covariance before real cross-study pooling.
- Add raw-data Nextflow pipelines for methylation, proteomics and metabolomics, following `workflows/rnaseq`.
- Fetch versioned references inside the RNA-seq pipeline instead of passing prepared files.
- Add R `metafor` cross-checks.
- Add richer forest plots.
- Add species-specific mapping releases.
- Add database export targets.
