#!/usr/bin/env bash
# DESeq2 for PRJNA913164 (QuantSeq 3' tag-seq): each heatwave treatment against the shared
# 20 C controls, blocking on ploidy. Salmon was run with --noLengthCorrection, and --tag-seq=yes
# drops the tximport transcript-length offset, since tag-seq yields one read per transcript.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 ANALYSIS_ROOT [REFERENCE_DIR]" >&2
  exit 2
fi

ANALYSIS_ROOT=$1
REFERENCE_DIR=${2:-"$ANALYSIS_ROOT/reference"}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
QUANT_DIR="$ANALYSIS_ROOT/salmon"
TX2GENE="$REFERENCE_DIR/GCF_963853765.1_tx2gene.tsv"

python "$SCRIPT_DIR/summarize_salmon_qc.py" --quant-dir "$QUANT_DIR" \
  --design-sheet "$ANALYSIS_ROOT/deseq2_samplesheet.csv" --output "$QUANT_DIR/salmon_qc_summary.tsv"
for test in heat heat_emersion; do
  case $test in
    heat) title="30 C seawater versus 20 C control (adult ctenidium)" ;;
    heat_emersion) title="30 C seawater then 6 h emersion at 44 C versus 20 C control (adult ctenidium)" ;;
  esac
  Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_$test.csv" \
    "$TX2GENE" "$ANALYSIS_ROOT/deseq2_$test" --test="$test" --replicates=0 --covariate=ploidy --tag-seq=yes \
    "--ma-title=$title"
done
