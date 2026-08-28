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
