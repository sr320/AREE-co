# Quantitative Proteomics Pipeline

Nextflow DSL2 workflow normalizing LFQ/TMT protein abundance tables, fitting empirical-Bayes moderated linear models (`limma::eBayes(robust=TRUE)`), and exporting AREE-validated protein evidence (`higher` / `lower`).

## Quick Test

```bash
python3 workflows/proteomics/tests/make_test_data.py --output /tmp/aree_prot_test
nextflow run workflows/proteomics/main.nf -profile test \
  --study_id SYNTH_PROT_TEST \
  --abundance_table /tmp/aree_prot_test/proteomics_abundance.tsv \
  --design /tmp/aree_prot_test/design.csv \
  --test_level heat --reference_level control \
  --outdir /tmp/aree_prot_out
```
