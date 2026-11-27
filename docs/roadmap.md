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

## Production Next Steps

- Expand the verified public-study collection, prioritizing measured resilience outcomes and independent replication.
- Add raw-data Nextflow pipelines for methylation, proteomics and metabolomics, following `workflows/rnaseq`.
- Fetch versioned references inside the RNA-seq pipeline instead of passing prepared files.
- Add richer forest plots.
- Add species-specific mapping releases.
- Add database export targets.
