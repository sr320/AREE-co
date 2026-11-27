# Metabolomics Differential Pipeline

Nextflow DSL2 workflow applying Probabilistic Quotient Normalization (PQN), `log2` transformation, empirical-Bayes moderated linear modeling (`limma`), and Wilcoxon rank-sum side checks, exporting AREE-validated metabolite evidence (`higher` / `lower`).

## Quick Test

```bash
python3 workflows/metabolomics/tests/make_test_data.py --output /tmp/aree_metab_test
nextflow run workflows/metabolomics/main.nf -profile test \
  --study_id SYNTH_METAB_TEST \
  --abundance_table /tmp/aree_metab_test/metabolomics_abundance.tsv \
  --design /tmp/aree_metab_test/design.csv \
  --test_level heat --reference_level control \
  --outdir /tmp/aree_metab_out
```
