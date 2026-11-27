# AREE Containers

Version-pinned container definitions for running Nextflow workflows and CLI tasks reproducibly with Docker or Singularity/Apptainer.

## Images

| Dockerfile | Default Tag | Contents | Used By |
| :--- | :--- | :--- | :--- |
| `Dockerfile.python` | `ghcr.io/sr320/aree-python:0.1.0` | Python 3.11 + `aree` CLI (`[site]` extras incl. `pyarrow`) | `FETCH_REFERENCE`, `EXPORT_EVIDENCE`, `VALIDATE_EVIDENCE` |
| `Dockerfile.rnaseq` | `ghcr.io/sr320/aree-rnaseq:0.1.0` | `fastqc=0.12.1`, `salmon=1.10.3`, R 4.3, `tximport=1.30.0`, `DESeq2=1.42.0` | `workflows/rnaseq` (`FASTQC`, `SALMON_*`, `DESEQ2`) |
| `Dockerfile.methylation` | `ghcr.io/sr320/aree-methylation:0.1.0` | `bismark`, `bowtie2`, `samtools`, R 4.3, `bioconductor-limma` | `workflows/methylation` (`BISMARK_*`, `LIMMA_DM`) |
| `Dockerfile.omics_stats` | `ghcr.io/sr320/aree-omics-stats:0.1.0` | R 4.3 (`limma`, `metafor`, `lme4`) + `aree` Python CLI | `workflows/proteomics`, `workflows/metabolomics` |

Every container includes `procps` (required by Nextflow process metrics and trace capture).

## Building Locally

Run from the repository root (`AREE-co/`):

```bash
docker build -f containers/Dockerfile.python -t ghcr.io/sr320/aree-python:0.1.0 .
docker build -f containers/Dockerfile.rnaseq -t ghcr.io/sr320/aree-rnaseq:0.1.0 .
docker build -f containers/Dockerfile.methylation -t ghcr.io/sr320/aree-methylation:0.1.0 .
docker build -f containers/Dockerfile.omics_stats -t ghcr.io/sr320/aree-omics-stats:0.1.0 .
```

## Running Nextflow with Containers

```bash
nextflow run workflows/rnaseq/main.nf -profile docker,test ...
nextflow run workflows/methylation/main.nf -profile singularity ...
```

Override image tags on the command line via `--container_python`, `--container_rnaseq`, `--container_methylation`, or `--container_omics_stats`.
