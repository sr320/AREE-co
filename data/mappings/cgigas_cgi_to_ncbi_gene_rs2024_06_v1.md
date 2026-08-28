# `cgigas_cgi_to_ncbi_gene_rs2024_06_v1`

## Scope and result

This mapping release covers the 150 significant features in study
`CGIG_HEAT_RNASEQ_PRJNA516762`. It maps legacy oyster_v9 `CGI_*` models to
NCBI Gene identifiers present in the *Crassostrea gigas* RefSeq annotation
`GCF_963853765.1-RS_2024_06`.

- 60 features have conservative, study-unique mappings with confidence
  `inferred`.
- 90 features remain `unresolved` and retain their original identifiers.
- All 60 resolved targets are unique within the selected feature set.
- No mapping is labeled `exact`, because assembly/annotation translation is
  involved even where the legacy coordinates match an NCBI release 101 model.

The machine-readable release is
`data/mappings/cgigas_cgi_to_ncbi_gene_rs2024_06_v1.tsv`.

## Method

1. Read the selected legacy feature IDs from the processed study table.
2. Recover legacy gene coordinates from the Ensembl Metazoa release 40
   oyster_v9 GFF3.
3. Translate oyster_v9 scaffold names to RefSeq sequence accessions using the
   `GCF_000297895.1` assembly report.
4. Select a mapping only when exactly one same-strand NCBI annotation release
   101 gene overlaps the legacy model and either legacy-model coverage is at
   least 0.8 or the smaller reciprocal coverage is at least 0.5.
5. Confirm that the NCBI Gene ID is present in the current RefSeq annotation.
   For discontinued Gene IDs, follow the NCBI Gene `currentid` field and then
   require the replacement to be present in the current annotation.
6. Downgrade ambiguity, absent models, missing replacements, and selected
   legacy-feature collisions onto one current Gene ID to `unresolved`.

The reproducible builder is `scripts/build_cgigas_identifier_map.py`.

## Evidence breakdown

| Outcome | Features |
| --- | ---: |
| Unique coordinate overlap; Gene ID remains current | 45 |
| Unique coordinate overlap; discontinued Gene ID has a current replacement | 15 |
| Legacy/new feature has no archived genomic model | 47 |
| No unique overlap meeting the coverage threshold | 24 |
| Discontinued Gene ID has no current replacement | 14 |
| Replacement Gene ID is absent from the current annotation | 2 |
| Multiple selected legacy features map to one current Gene ID | 2 |
| Multiple qualifying legacy-annotation overlaps | 1 |

## Source releases and checksums

- Ensembl Metazoa release 40 oyster_v9 GFF3:
  <https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-40/gff3/crassostrea_gigas/Crassostrea_gigas.oyster_v9.40.gff3.gz>
  - SHA-256 of decompressed GFF3:
    `c435c0bd458f47680c51de87f4eaf279b91fde9c8cd312101b89dd23f6dd9adb`
  - Published MD5 of decompressed GFF3:
    `90a747fbc94a0a9225c43f75cc40b9db`
- NCBI oyster_v9 assembly `GCF_000297895.1` and annotation release 101:
  <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000297895.1/>
  - Assembly report SHA-256:
    `79d8f326dde488d721bf2487e464fa6c27eb6081d8068d8991ae99abb9c032d8`
  - Genomic GFF3 SHA-256:
    `321cd0c3e2c784b9c6ed208afe61eb8a05d236eb5a4b0971edd4388c9f066b65`
- NCBI RefSeq assembly `GCF_963853765.1`, build xbMagGiga1.1,
  annotation `RS_2024_06`:
  <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_963853765.1/>
  - Genomic GFF3 SHA-256:
    `18990f45004e2db2c607be33f1d3f7512bb721d90fd272f489f6c9e64ee0b92e`
- NCBI Gene E-utilities summaries for discontinued candidate Gene IDs:
  <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
  - Retrieved summary JSON SHA-256:
    `ac21e74c889e0a74915e2b4b6862ace76864961250ca37181a5a636666b21997`

Source annotations were downloaded on 2026-08-27. The large source files are
not committed; checksums and the builder make the mapping release auditable.

## Limitations

This is an annotation-translation release, not a sequence-alignment liftover.
Coordinate overlap is evaluated only on the shared oyster_v9 assembly before
current NCBI Gene status is checked. Unresolved features must not be collapsed
across studies without additional sequence- or transcript-level evidence.
