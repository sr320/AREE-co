# Reanalysis intake: PRJNA694496

## Why this study is the second analytical priority

Tan et al. (2021) compare constitutive larval RNA expression between a Pacific
oyster population selected for thermotolerance and a non-selected control.
Unlike the first AREE real-study intake, this study links the molecular contrast
to a favorable survival phenotype: 78.05% versus 66.34% after acute heat shock
and 73% versus 62% during the field summer trial.

The study has six public paired-end RNA-seq runs (three per population), but the
supplement contains pathway summaries rather than the complete 4,764-gene
differential-expression table. AREE therefore treats it as a raw-reanalysis
target, not as processed evidence ready for harmonization.

## Verified run assignment

NCBI stores all six runs under one BioSample and does not encode the selected
and control labels directly. Supplementary Figure S2 reports a distinct clean
read count for each named replicate. Those six counts exactly match the ENA run
records, resolving the assignment without guessing.

| Condition | Replicate | SRA run | SRA library | Clean read pairs |
| --- | --- | --- | --- | ---: |
| Selected | A | `SRR13576947` | `N11A1` | 23,941,716 |
| Selected | B | `SRR13576945` | `N11B1` | 20,502,956 |
| Selected | C | `SRR13576944` | `N11C1` | 22,183,461 |
| Control | A | `SRR13576943` | `N4A1` | 21,776,617 |
| Control | B | `SRR13576942` | `N4B1` | 21,642,163 |
| Control | C | `SRR13576946` | `N4C1` | 22,367,076 |
| **Total** |  |  |  | **132,413,989** |

The tracked run manifest includes paired FASTQ URLs, byte counts, and MD5
checksums from ENA.

## Reanalysis contract

- Reference: *C. gigas* RefSeq assembly `GCF_963853765.1`, annotation
  `RS_2024_06`, matching AREE's current identifier target.
- Contrast: selected minus control. Positive log2 fold change means higher
  constitutive expression in the thermotolerance-selected population.
- Replication: three RNA-seq libraries per population.
- Quantification: transcript-aware quantification followed by gene-level
  summarization against the versioned annotation.
- Differential expression: DESeq2 with all tested genes exported, including
  base mean, log2 fold change, `lfcSE`, Wald statistic, raw p-value, and adjusted
  p-value.
- Sensitivity: report both unshrunk coefficients (for standard errors and
  meta-analysis) and shrunken log2 fold changes (for ranking only).
- QC gates: verify FASTQ MD5, library layout, read yield, mapping/assignment
  rate, sample correlation, PCA separation, and outlier diagnostics before
  evidence harmonization.
- Identifier handling: emit current NCBI Gene IDs directly and preserve the
  source transcript IDs; do not reuse the legacy coordinate mapping for
  current-reference output.

## Resource preflight

The paired FASTQs total 16,252,181,463 compressed bytes (16.25 GB decimal).
The repository filesystem is too small for the inputs, so raw data and workflow
artifacts are staged under `/Volumes/omics/scratch/AREE/`. The storage preflight
found 490.01 GB available there against a 30.07 GB direct-from-gzip requirement.
Retain at least 60 GB if trimmed FASTQ copies are enabled later.

The resumable, MD5-verified staging command is:

```bash
python scripts/prepare_prjna694496_fastqs.py \
  --output-dir /Volumes/omics/scratch/AREE/CGIG_THERMOTOL_RNASEQ_PRJNA694496/fastq \
  --write-samplesheet /Volumes/omics/scratch/AREE/CGIG_THERMOTOL_RNASEQ_PRJNA694496/nfcore_samplesheet.csv \
  --write-design-sheet /Volumes/omics/scratch/AREE/CGIG_THERMOTOL_RNASEQ_PRJNA694496/deseq2_samplesheet.csv \
  --workers 12 \
  --download
```

The script fetches and MD5-verifies the six pairs and creates separate input
sheets for quantification and DESeq2 design metadata. The originally pinned
nf-core/rnaseq 3.26.0 target cannot run reproducibly on this host: its installed
Nextflow 20.10.0 is obsolete, no container runtime is available, and nf-core
does not support its Conda profile on macOS. The local execution path therefore
uses a versioned native Salmon/tximport/DESeq2 workflow instead of silently
relaxing those requirements.

Stage and verify the current RefSeq genome, transcriptome, annotation, and
transcript-to-GeneID map with:

```bash
python scripts/prepare_gcf963853765_reference.py \
  --output-dir /Volumes/omics/scratch/AREE/CGIG_THERMOTOL_RNASEQ_PRJNA694496/reference
```

The tracked environment specification is
`config/CGIG_THERMOTOL_RNASEQ_PRJNA694496_analysis_environment.yaml`. The active
scratch environment pins Salmon 1.10.3, FastQC 0.12.1, tximport 1.30.0, DESeq2
1.42.0, and apeglm 1.24.0. Run the complete QC, Salmon quantification, tximport,
and selected-versus-control DESeq2 analysis with:

```bash
micromamba run \
  --prefix /Volumes/omics/scratch/AREE/CGIG_THERMOTOL_RNASEQ_PRJNA694496/analysis-env \
  bash scripts/run_prjna694496_salmon.sh \
  /Volumes/omics/scratch/AREE/CGIG_THERMOTOL_RNASEQ_PRJNA694496 8
```

The workflow exports all tested gene-level effects and standard errors, an
apeglm ranking table, normalized counts, PCA coordinates and plot, sample
correlations, the serialized DESeq2 object, and session information.

## Interpretation guardrail

This is a population-selection association. RNA was collected from larval
pools at ambient conditions, while survival was assessed in related later-stage
cohorts. The study can support repeated resilience-associated genes, but it
does not establish that larval expression predicts individual survival.

## Completion condition

The study becomes meta-analysis eligible only after raw reanalysis produces
versioned gene-level effects and standard errors, QC passes, and at least one
current Gene ID overlaps the first real-study release. Because the contrasts
differ (constitutive selected-versus-control versus acute heat response), any
pooled model must stratify by contrast class or explicitly model that
moderator.
