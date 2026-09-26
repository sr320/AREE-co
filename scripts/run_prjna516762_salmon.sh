#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 4 ]]; then
  echo "usage: $0 ANALYSIS_ROOT [THREADS] [SAMPLESHEET] [REFERENCE_DIR]" >&2
  exit 2
fi

ANALYSIS_ROOT=$1
THREADS=${2:-8}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SAMPLESHEET=${3:-"$ANALYSIS_ROOT/nfcore_samplesheet.csv"}
REFERENCE_DIR=${4:-"$ANALYSIS_ROOT/reference"}
FASTQ_DIR="$ANALYSIS_ROOT/fastq"
QUANT_DIR="$ANALYSIS_ROOT/salmon"
QC_DIR="$ANALYSIS_ROOT/fastqc"
RESULTS_DIR="$ANALYSIS_ROOT/deseq2"
DESIGN_SHEET="$ANALYSIS_ROOT/deseq2_samplesheet.csv"
TX2GENE="$REFERENCE_DIR/GCF_963853765.1_tx2gene.tsv"
SALMON_INDEX="$REFERENCE_DIR/salmon_index_decoy"

for executable in salmon fastqc Rscript; do
  command -v "$executable" >/dev/null 2>&1 || { echo "missing required executable: $executable" >&2; exit 1; }
done
for input_path in "$SAMPLESHEET" "$DESIGN_SHEET" "$TX2GENE" "$SALMON_INDEX/versionInfo.json"; do
  [[ -s "$input_path" ]] || { echo "missing required input: $input_path" >&2; exit 1; }
done

mkdir -p "$QUANT_DIR" "$QC_DIR" "$RESULTS_DIR"
missing_fastqc=()
while IFS=, read -r sample fastq_1 fastq_2 strandedness; do
  for fastq_path in "$fastq_1" "$fastq_2"; do
    [[ -s "$fastq_path" ]] || { echo "missing FASTQ: $fastq_path" >&2; exit 1; }
    fastq_name=$(basename "$fastq_path" .fastq.gz)
    [[ -s "$QC_DIR/${fastq_name}_fastqc.zip" ]] || missing_fastqc+=("$fastq_path")
  done
done < <(tail -n +2 "$SAMPLESHEET")
if [[ ${#missing_fastqc[@]} -gt 0 ]]; then
  fastqc --threads "$THREADS" --outdir "$QC_DIR" "${missing_fastqc[@]}"
fi

tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness; do
  sample_quant="$QUANT_DIR/$sample"
  if [[ -s "$sample_quant/quant.sf" && -s "$sample_quant/aux_info/meta_info.json" ]] &&
    python -c 'import json,sys; m=json.load(open(sys.argv[1])); sys.exit(0 if m.get("end_time") and not m.get("quant_errors") and m.get("num_processed",0)>0 else 1)' "$sample_quant/aux_info/meta_info.json"; then
    echo "verified existing Salmon result: $sample"
    continue
  fi
  salmon quant --index "$SALMON_INDEX" --libType A --mates1 "$fastq_1" --mates2 "$fastq_2" \
    --threads "$THREADS" --validateMappings --seqBias --gcBias --output "$sample_quant"
done

python "$SCRIPT_DIR/summarize_salmon_qc.py" --quant-dir "$QUANT_DIR" \
  --design-sheet "$DESIGN_SHEET" --output "$QUANT_DIR/salmon_qc_summary.tsv"
Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$DESIGN_SHEET" "$TX2GENE" "$RESULTS_DIR" \
  --test=heat "--ma-title=Heat shock versus control"
