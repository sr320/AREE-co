# Identifier Hierarchy

AREE preserves original feature identifiers and adds standardized identifiers using this hierarchy:

1. Curated NCBI Gene or Ensembl gene ID from the study species and genome version.
2. One-to-one ortholog reference gene for cross-species comparison.
3. UniProt accession linked to a gene model.
4. Transcript, protein, or locus ID linked to a gene model.
5. Orthogroup or inferred feature class.
6. Original ID retained with `unresolved` mapping confidence.

Mapping confidence values are: `exact`, `one-to-one ortholog`, `one-to-many ortholog`, `many-to-one ortholog`, `inferred`, and `unresolved`.

