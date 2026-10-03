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
  # A single-end sample leaves fastq_2 empty.
  for fastq_path in "$fastq_1" ${fastq_2:+"$fastq_2"}; do
    [[ -s "$fastq_path" ]] || { echo "missing FASTQ: $fastq_path" >&2; exit 1; }
    fastq_name=$(basename "$fastq_path" .fastq.gz)
    [[ -s "$QC_DIR/${fastq_name}_fastqc.zip" ]] || missing_fastqc+=("$fastq_path")
  done
done < <(tail -n +2 "$SAMPLESHEET")
# Nothing downstream reads the FastQC archives -- the recorded quality_control_status is
# built from Salmon mapping rates and the DESeq2 summary -- so on a study large enough for
# FastQC to take hours, SKIP_FASTQC=1 defers it and lets the quantification run first.
if [[ ${#missing_fastqc[@]} -gt 0 && "${SKIP_FASTQC:-0}" != "1" ]]; then
  fastqc --threads "$THREADS" --outdir "$QC_DIR" "${missing_fastqc[@]}"
elif [[ ${#missing_fastqc[@]} -gt 0 ]]; then
  echo "SKIP_FASTQC=1: deferring FastQC for ${#missing_fastqc[@]} files"
fi

tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq_1 fastq_2 strandedness; do
  sample_quant="$QUANT_DIR/$sample"
  if [[ -s "$sample_quant/quant.sf" && -s "$sample_quant/aux_info/meta_info.json" ]] &&
    python -c 'import json,sys; m=json.load(open(sys.argv[1])); sys.exit(0 if m.get("end_time") and not m.get("quant_errors") and m.get("num_processed",0)>0 else 1)' "$sample_quant/aux_info/meta_info.json"; then
    echo "verified existing Salmon result: $sample"
    continue
  fi
  if [[ -n "$fastq_2" ]]; then
    reads=(--mates1 "$fastq_1" --mates2 "$fastq_2" --gcBias)
  else
    # Salmon's GC-bias model is built on fragment lengths, which single-end reads do not
    # observe, so single-end samples are corrected for sequence bias only.
    reads=(--unmatedReads "$fastq_1")
  fi
  salmon quant --index "$SALMON_INDEX" --libType A "${reads[@]}" \
    --threads "$THREADS" --validateMappings --seqBias --output "$sample_quant"
done

python "$SCRIPT_DIR/summarize_salmon_qc.py" --quant-dir "$QUANT_DIR" \
  --design-sheet "$DESIGN_SHEET" --output "$QUANT_DIR/salmon_qc_summary.tsv"
# Six pools per temperature, three in each of the two full-sib families, all sampled at 0 h. Each
# pool is 10 different spat from one tank, as in the starvation arm whose pools passed the
# expressed-SNP independence check, so they are not collapsed; family is blocked out as
# ~ family + condition. The 23 C pools are the starvation arm's fed controls.
Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$DESIGN_SHEET" "$TX2GENE" "$RESULTS_DIR" \
  --test=heat --replicates=6 --covariate=family \
  "--ma-title=30 C versus 23 C, fed (0 h)"
