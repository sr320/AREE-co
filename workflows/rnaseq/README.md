# AREE RNA-seq reanalysis pipeline

Paired-end FASTQs → FastQC → Salmon (decoy-aware selective alignment) → tximport/DESeq2 → AREE processed
evidence, exact identifier mappings and RefSeq annotations, ready for `aree harmonize`.

Reads are not trimmed; Salmon's selective alignment soft-clips adapters and low-quality ends, matching the
AREE reanalyses of PRJNA694496 and PRJNA516762. Tested with Nextflow 26.04 using both Salmon 1.10.3 and
2.3.4, DESeq2 1.42 and tximport 1.30.

## Requirements

- Nextflow ≥ 24.04 and Java 17+.
- `aree` installed in the environment that launches Nextflow (`pip install -e ".[dev]"` from the repo root).
- Either `salmon`, `fastqc` and `Rscript` with DESeq2/tximport on `PATH` (default profile), or
  micromamba to build them from `environment.yml` (`-profile micromamba`).

## Inputs

| Parameter | Content |
|---|---|
| `--study_id` | Registry study ID; names the evidence files. |
| `--samplesheet` | CSV `sample,fastq_1,fastq_2[,strandedness]` (written by `scripts/prepare_*_fastqs.py`). |
| `--design` | CSV `sample,condition,replicate,run_accession`; `sample` must match the sample sheet. |
| `--tx2gene` | TSV `transcript_id,gene_id`, gene IDs as `NCBI:GeneID:<id>`. |
| `--gff` | Gzipped RefSeq GFF of the same release, for gene annotations. |
| `--transcripts` [+ `--decoys`] or `--salmon_index` | Gentrome and decoy list to index, or a prebuilt index. |
| `--test_level` | Design condition compared against `--reference_level` (default `control`). |

Optional: `--replicates` (exact count per condition; default 0 = any ≥ 2), `--min_samples` (default 3),
`--quality_flags` (comma-separated study flags), `--sample_comparison`, `--reference_release`,
`--skip_fastqc`, `--outdir` (default `results/rnaseq/<study_id>`).

## Example: PRJNA516762 heat shock

```bash
python scripts/prepare_gcf963853765_reference.py            # GCF_963853765.1 gentrome, decoys, tx2gene
python scripts/prepare_prjna516762_fastqs.py --analysis-root results/PRJNA516762 --download
nextflow run workflows/rnaseq -profile micromamba \
  --study_id CGIG_HEAT_RNASEQ_PRJNA516762 --test_level heat --replicates 3 \
  --samplesheet results/PRJNA516762/nfcore_samplesheet.csv \
  --design results/PRJNA516762/deseq2_samplesheet.csv \
  --transcripts data/reference/GCF_963853765.1/GCF_963853765.1_gentrome.fa \
  --decoys data/reference/GCF_963853765.1/GCF_963853765.1_decoys.txt \
  --tx2gene data/reference/GCF_963853765.1/GCF_963853765.1_tx2gene.tsv \
  --gff data/reference/GCF_963853765.1/GCF_963853765.1_xbMagGiga1.1_genomic.gff.gz \
  --quality_flags raw_reanalysis,acute_heat_response,pooled_libraries \
  --sample_comparison heat_35C_2h_vs_control_12C \
  --ma_title "Heat shock versus control"
aree harmonize --study CGIG_HEAT_RNASEQ_PRJNA516762 \
  --input results/rnaseq/CGIG_HEAT_RNASEQ_PRJNA516762/evidence/CGIG_HEAT_RNASEQ_PRJNA516762_rnaseq.tsv \
  --mapping results/rnaseq/CGIG_HEAT_RNASEQ_PRJNA516762/evidence/CGIG_HEAT_RNASEQ_PRJNA516762_mapping.tsv \
  --output data/harmonized/evidence.tsv
```

The per-study shell wrappers in `scripts/run_prjna*_salmon.sh` record how the committed reanalyses were run;
new studies should use this pipeline.

## Outputs (`--outdir`)

- `fastqc/`: per-FASTQ reports.
- `salmon/<sample>/`: Salmon quantifications; `salmon/salmon_qc_summary.tsv`: mapping rates in design order.
- `deseq2/`: `<test>_vs_<reference>_deseq2_all_genes.tsv` (unshrunk effects used as evidence), normalized
  counts, PCA, VST correlations, Cook's-distance diagnostics, MA and PCA plots, `sessionInfo.txt`.
- `evidence/`: `<study_id>_rnaseq.tsv`, `<study_id>_mapping.tsv`, `<study_id>_annotations.tsv`.

## Testing

`workflows/rnaseq/tests/make_test_data.py` writes a deterministic synthetic study (300 genes, 3 control and
3 heat libraries; genes 0–29 simulated 4× up and 30–59 4× down). `tests/test_nextflow_rnaseq.py` runs the
pipeline on it, with `-stub-run` and for real, and checks that DESeq2 recovers the simulated changes and that
the evidence harmonizes. It skips unless the tools are available; set `AREE_NEXTFLOW` to choose a Nextflow
launcher and `AREE_REQUIRE_PIPELINE=1` to fail instead of skip.
