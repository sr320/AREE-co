# WGBS / Methylation Differential Pipeline

Nextflow DSL2 workflow fitting empirical-Bayes M-value linear models (`limma`) on Bismark CpG coverage (`_M`/`_U` counts or beta matrices) and exporting AREE-validated methylation evidence (`hypermethylated` / `hypomethylated`).

## Quick Test

```bash
python3 workflows/methylation/tests/make_test_data.py --output /tmp/aree_meth_test
nextflow run workflows/methylation/main.nf -profile test \
  --study_id SYNTH_METH_TEST \
  --coverage_table /tmp/aree_meth_test/methylation_coverage.tsv \
  --design /tmp/aree_meth_test/design.csv \
  --test_level heat --reference_level control \
  --outdir /tmp/aree_meth_out
```
