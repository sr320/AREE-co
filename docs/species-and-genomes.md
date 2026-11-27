# Adding Species and Handling Genome Versions

To add a species:

1. Add species-specific genome and annotation metadata to study YAML files.
2. Add identifier maps under `data/mappings/`.
3. Preserve original locus, transcript, protein, and metabolite identifiers.
4. Add ortholog links only when confidence can be documented.
5. Record mapping confidence for every mapped feature.

Genome-version changes should be treated as provenance events. Do not overwrite historical standardized identifiers without recording the mapping release.

## Pacific oyster reference used by the first real-study intake

Mapping release `cgigas_cgi_to_ncbi_gene_rs2024_06_v1` translates selected
oyster_v9 (`GCA_000297895.1`) features to NCBI Gene IDs confirmed in RefSeq
assembly `GCF_963853765.1`, annotation `RS_2024_06`. The release table and its
method, source checksums, thresholds, and unresolved cases are stored together
under `data/mappings/`.

## Versioned Multi-Species Bivalve Reference Catalog & Ortholog Releases

Eight RefSeq bivalve assemblies are cataloged in `config/references.yaml` and fetched/verified via `aree fetch-reference` (and Nextflow module `modules/fetch_reference.nf`):

| Code | Species | Assembly | Annotation Release | TaxID |
| :--- | :--- | :--- | :--- | :--- |
| `mgigas` | *Magallana gigas* (*Crassostrea gigas*) | `GCF_963853765.1` | `GCF_963853765.1-RS_2024_06` | 29159 |
| `cvirg` | *Crassostrea virginica* | `GCF_053477285.1` | `GCF_053477285.1-RS_2025_12` | 6565 |
| `medul` | *Mytilus edulis* | `GCF_963676685.1` | `GCF_963676685.1-RS_2025_02` | 6550 |
| `mgall` | *Mytilus galloprovincialis* | `GCF_965363235.1` | `GCF_965363235.1-RS_2025_08` | 29158 |
| `mcali` | *Mytilus californianus* | `GCF_021869535.1` | `NCBI Mytilus californianus Annotation Release 100` | 6549 |
| `mtros` | *Mytilus trossulus* | `GCF_036588685.1` | `GCF_036588685.1-RS_2024_02` | 6551 |
| `mmerc` | *Mercenaria mercenaria* | `GCF_021730395.1` | `GCF_021730395.1-RS_2023_01` | 6596 |
| `rphil` | *Ruditapes philippinarum* | `GCF_026571515.1` | `GCF_026571515.1-RS_2023_10` | 129788 |

OrthoFinder 3.1.5 pairwise ortholog maps (`orthofinder3.1.5_8bivalves_2026-11-26_v1`) linking all seven partner species to *Magallana gigas* GeneIDs live under `data/mappings/orthologs/` (rebuilt by `scripts/build_multispecies_ortholog_maps.py`).
