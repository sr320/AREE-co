# Adding Species and Handling Genome Versions

To add a species:

1. Add species-specific genome and annotation metadata to study YAML files.
2. Add identifier maps under `data/mappings/`.
3. Preserve original locus, transcript, protein, and metabolite identifiers.
4. Add ortholog links only when confidence can be documented.
5. Record mapping confidence for every mapped feature.

Genome-version changes should be treated as provenance events. Do not overwrite historical standardized identifiers without recording the mapping release.

