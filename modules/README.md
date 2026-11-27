# Nextflow Modules

Shared DSL2 process modules reused across AREE-co workflows (`workflows/rnaseq/`, `workflows/methylation/`, `workflows/proteomics/`, `workflows/metabolomics/`):

- `fetch_reference.nf` (`FETCH_REFERENCE`): Fetches and verifies a versioned multi-species RefSeq bundle from `config/references.yaml` via `aree fetch-reference`, generating `tx2gene.tsv`, `genomic.gff.gz`, and optional decoy-aware `gentrome.fa.gz` + `decoys.txt`.
- `validate_evidence.nf` (`VALIDATE_EVIDENCE`): Validates exported `processed_evidence.tsv`, `identifier_mapping.tsv`, and `feature_annotations.tsv` against AREE schema constraints via `aree validate-processed`.
