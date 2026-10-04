# Adding a Study

1. Create `registry/studies/STUDY_ID.yaml` from `registry/study_template.yaml`.
2. Use controlled vocabulary terms for `phenotype` and `stressor_class`.
3. Preserve original treatment and control labels.
4. Record raw and processed data availability separately.
5. Add caveats instead of hiding missing metadata.
6. Validate the file:

```bash
aree validate-study registry/studies/STUDY_ID.yaml
```

7. Add it to the registry:

```bash
aree register-study registry/studies/STUDY_ID.yaml
```

## Raw Reanalysis Versus Processed Harmonization

Set `sample_size` to the count included in the final primary analysis, after sample exclusions or technical-replicate collapse. Use `biological_replication` to identify the counted unit (individual animals, pools or libraries), counts per group, exclusions, and the independent experimental units. For example, PRJEB18614 uses 24 pools after excluding two of 26 libraries; PRJNA735889 uses 50 animals after excluding the resequenced 13Ebis library. A pooled-library count is not an individual-animal or independent-tank count.

After changing the primary analysis or its registry count, rerun harmonization with the matching mapping release before scoring. The committed-evidence consistency tests compare every row's count with both YAML and CSV metadata, and scoring rejects conflicting counts within a study rather than choosing one silently.

Use raw-data reanalysis when public FASTQ, spectra, or feature-level raw files are available and licensing permits reuse. Use processed-results harmonization when only publication supplements or repository-derived result tables are available.

Write real-study evidence outside the synthetic demo directory:

```bash
aree harmonize --study STUDY_ID --input data/processed/STUDY_ID_rnaseq.tsv --mapping data/mappings/MAPPING_RELEASE.tsv --output data/harmonized
```

`--output` is required for any study whose `data_availability.status` is not `simulated_*`, so real evidence cannot land in the demo table by accident. Given a directory, `harmonize` writes one `STUDY_ID.tsv` per study (rerunning a study replaces only its own file), which keeps every committed file well under GitHub's 50 MB warning and 100 MB limit. Point the downstream commands at the same directory; they concatenate every `*.tsv` in it (a single TSV also works). Without `--evidence` they read the demo data:

```bash
aree meta-analyze --evidence data/harmonized --output results/meta_analysis.tsv
aree score-candidates --evidence data/harmonized --output results/candidate_scores.tsv
aree build-evidence-cards --evidence data/harmonized --scores results/candidate_scores.tsv --output-dir results/evidence_cards
```

Scoring recomputes heterogeneity from the evidence it scores, so `meta_analysis.tsv` is for inspection only and a filtered or stale copy cannot change the ranking. Effects are pooled and summarized only within one `effect_size_type`; candidates remain one row per feature so cross-context and multi-omics convergence still count. Use `--phenotype`/`--stressor` on `score-candidates` for a context-specific ranking.

Rerunning `aree harmonize` on unchanged inputs leaves the evidence table byte-identical (`date_generated` is kept); it only changes when the records do.
