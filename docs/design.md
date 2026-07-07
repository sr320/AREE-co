# Phase 1 Design

## Assumptions

- Public studies will vary in raw-data availability, so AREE supports raw-data reanalysis and processed-results harmonization.
- The first production species is *Crassostrea gigas*, but evidence records include species and ortholog fields from the start.
- Identifier mappings are versioned data products with confidence levels, not silent replacements.
- Phenotype labels need controlled terms while preserving original study descriptions.
- Candidate ranking must be transparent and tunable, not a black-box model.

## Architecture

AREE has five connected layers:

1. Study registry and dataset intake.
2. Standardized reanalysis workflow scaffolds.
3. Cross-study evidence harmonization.
4. Meta-analysis and candidate prioritization.
5. User-facing reports and interface.

## Provenance Model

Each evidence row retains study ID, source accession, input file, checksum, workflow version, software method, genome and annotation version, generation date, quality flags, and curation caveats.

## Ranking Framework

The MVP candidate score combines study count, biological sample size, effect magnitude, adjusted significance, direction consistency, phenotype relevance, tissue/life-stage breadth, assay diversity, mapping confidence, quality flags, and heterogeneity penalty.

