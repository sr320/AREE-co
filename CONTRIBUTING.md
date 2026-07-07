# Contributing

AREE welcomes code, documentation, workflow modules, and curated public-study metadata.

## Add a Study

1. Copy `registry/studies/CGIG_HEAT_RNASEQ_001.yaml`.
2. Replace identifiers, accessions, study context, assay metadata, and caveats.
3. Preserve original treatment labels and map them separately to controlled terms.
4. Run `aree validate-study registry/studies/YOUR_STUDY.yaml`.
5. Run `aree register-study registry/studies/YOUR_STUDY.yaml`.

## Evidence Claims

AREE reports associations and evidence convergence. Do not label a candidate as validated unless validation evidence exists outside this resource and is cited.

## Data Governance

Only register public or appropriately shareable datasets. Preserve provenance links, raw/processed availability, limitations, and manual curation decisions.

