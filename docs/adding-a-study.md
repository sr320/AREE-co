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

Use raw-data reanalysis when public FASTQ, spectra, or feature-level raw files are available and licensing permits reuse. Use processed-results harmonization when only publication supplements or repository-derived result tables are available.

Write real-study evidence outside the synthetic demo directory:

```bash
aree harmonize --study STUDY_ID --input data/processed/STUDY_ID_rnaseq.tsv --mapping data/mappings/MAPPING_RELEASE.tsv --output data/harmonized/evidence.tsv
```

`--output` is required for any study whose `data_availability.status` is not `simulated_*`, so real evidence cannot land in the demo table by accident. Point the downstream commands at the same table; without `--evidence` they read the demo data:

```bash
aree meta-analyze --evidence data/harmonized/evidence.tsv --output results/meta_analysis.tsv
aree score-candidates --evidence data/harmonized/evidence.tsv --meta results/meta_analysis.tsv --output results/candidate_scores.tsv
aree build-evidence-cards --evidence data/harmonized/evidence.tsv --meta results/meta_analysis.tsv --scores results/candidate_scores.tsv --output-dir results/evidence_cards
```

Rerunning `aree harmonize` on unchanged inputs leaves the evidence table byte-identical (`date_generated` is kept); it only changes when the records do.
