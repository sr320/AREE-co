# Multi-Species Bivalve Ortholog Mapping Releases

Versioned gene-level ortholog mappings translating 7 aquaculture bivalve RefSeq annotations to *Magallana gigas* (`GCF_963853765.1`, annotation `GCF_963853765.1-RS_2024_06`), generated from OrthoFinder 3.1.5 (`orthofinder3.1.5_8bivalves_2026-11-26_v1`) across longest-isoform RefSeq primary proteomes.

## Files

- `mgigas_orthogroups_v1.tsv`: *Magallana gigas* GeneID-anchored orthogroup table (`25,290` indexed *M. gigas* genes across the 8 species).
- `<sp>_to_mgigas_orthofinder_v1.tsv`: 1-row-per-source-gene mapping tables compatible with `aree.harmonize.identifiers.load_mapping`:
  - `cvirg_to_mgigas_orthofinder_v1.tsv` (*Crassostrea virginica*, `GCF_053477285.1`, `22,213` genes)
  - `medul_to_mgigas_orthofinder_v1.tsv` (*Mytilus edulis*, `GCF_963676685.1`, `19,504` genes)
  - `mgall_to_mgigas_orthofinder_v1.tsv` (*Mytilus galloprovincialis*, `GCF_965363235.1`, `17,339` genes)
  - `mcali_to_mgigas_orthofinder_v1.tsv` (*Mytilus californianus*, `GCF_021869535.1`, `16,488` genes)
  - `mtros_to_mgigas_orthofinder_v1.tsv` (*Mytilus trossulus*, `GCF_036588685.1`, `19,154` genes)
  - `mmerc_to_mgigas_orthofinder_v1.tsv` (*Mercenaria mercenaria*, `GCF_021730395.1`, `18,872` genes)
  - `rphil_to_mgigas_orthofinder_v1.tsv` (*Ruditapes philippinarum*, `GCF_026571515.1`, `20,680` genes)
- `release_summary.tsv`: per-species mapping counts by cardinality (`one-to-one ortholog`, `many-to-one ortholog`, `one-to-many ortholog`, `inferred` for M:N).

## Cardinality Rules

Every `feature_id_original` appears at most once in each pairwise release table:
- `1:1` (`one-to-one ortholog`): single source gene and single *M. gigas* target gene.
- `M:1` (`many-to-one ortholog`): multiple source co-orthologs mapping to one *M. gigas* gene.
- `1:M` (`one-to-many ortholog`): one source gene mapping to multiple *M. gigas* co-orthologs; `feature_id_standardized` takes the lowest numeric *M. gigas* GeneID and `candidate_mgigas_gene_ids` preserves all targets.
- `M:N` (`inferred`): multiple source genes and multiple *M. gigas* co-orthologs; tagged `inferred` with all candidate *M. gigas* GeneIDs recorded.
