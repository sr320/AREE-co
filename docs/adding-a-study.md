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
aree harmonize --study STUDY_ID --input data/processed/STUDY_ID_rnaseq.tsv --output data/harmonized/evidence.tsv
```
