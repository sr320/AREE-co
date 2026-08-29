#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 ANALYSIS_ROOT [THREADS]" >&2
  exit 2
fi

ANALYSIS_ROOT=$1
THREADS=${2:-8}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ASSEMBLY=GCF_963853765.1_xbMagGiga1.1
REFERENCE_DIR="$ANALYSIS_ROOT/reference"
FASTQ_DIR="$ANALYSIS_ROOT/fastq"
QUANT_DIR="$ANALYSIS_ROOT/salmon"
QC_DIR="$ANALYSIS_ROOT/fastqc"
RESULTS_DIR="$ANALYSIS_ROOT/deseq2"
SAMPLESHEET="$ANALYSIS_ROOT/nfcore_samplesheet.csv"
DESIGN_SHEET="$ANALYSIS_ROOT/deseq2_samplesheet.csv"
TX2GENE="$REFERENCE_DIR/GCF_963853765.1_tx2gene.tsv"
GENTROME="$REFERENCE_DIR/GCF_963853765.1_gentrome.fa"
DECOYS="$REFERENCE_DIR/GCF_963853765.1_decoys.txt"
SALMON_INDEX="$REFERENCE_DIR/salmon_index_decoy"

for executable in salmon fastqc Rscript; do
  if ! command -v "$executable" >/dev/null 2>&1; then
    echo "missing required executable: $executable" >&2
    exit 1
  fi
done
for input_path in "$SAMPLESHEET" "$DESIGN_SHEET" "$GENTROME" "$DECOYS" "$TX2GENE"; do
  if [[ ! -s "$input_path" ]]; then
    echo "missing required input: $input_path" >&2
    exit 1
  fi
done

mkdir -p "$QUANT_DIR" "$QC_DIR" "$RESULTS_DIR"
if [[ ! -s "$SALMON_INDEX/versionInfo.json" ]]; then
  salmon index \
    --transcripts "$GENTROME" \
    --decoys "$DECOYS" \
    --keepDuplicates \
    --index "$SALMON_INDEX" \
    --threads "$THREADS"
fi

missing_fastqc=()
for fastq_path in "$FASTQ_DIR"/*.fastq.gz; do
  fastq_name=$(basename "$fastq_path" .fastq.gz)
  if [[ ! -s "$QC_DIR/${fastq_name}_fastqc.zip" ]]; then
    missing_fastqc+=("$fastq_path")
  fi
done
if [[ ${#missing_fastqc[@]} -gt 0 ]]; then
  fastqc --threads "$THREADS" --outdir "$QC_DIR" "${missing_fastqc[@]}"
fi

tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness; do
  sample_quant="$QUANT_DIR/$sample"
  if [[ -s "$sample_quant/quant.sf" ]]; then
    echo "verified existing Salmon result: $sample"
    continue
  fi
  salmon quant \
    --index "$SALMON_INDEX" \
    --libType A \
    --mates1 "$fastq_1" \
    --mates2 "$fastq_2" \
    --threads "$THREADS" \
    --validateMappings \
    --seqBias \
    --gcBias \
    --output "$sample_quant"
done

python "$SCRIPT_DIR/summarize_salmon_qc.py" \
  --quant-dir "$QUANT_DIR" \
  --design-sheet "$DESIGN_SHEET" \
  --output "$QUANT_DIR/salmon_qc_summary.tsv"

Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" \
  "$QUANT_DIR" \
  "$DESIGN_SHEET" \
  "$TX2GENE" \
  "$RESULTS_DIR"
