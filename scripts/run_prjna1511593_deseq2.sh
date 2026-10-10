#!/usr/bin/env bash
# DESeq2 for PRJNA1511593 (pure C. gigas gill): 37 C for 12 h and 42 C for 1 h, each against the
# shared controls, 3 v 3. The treatments are read from the sample names (H37H12, H42H1, H0).
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
for test in heat_37C_12h heat_42C_1h; do
  case $test in
    heat_37C_12h) title="37 C for 12 h versus control (gill)" ;;
    heat_42C_1h) title="42 C for 1 h versus control (gill)" ;;
  esac
  Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_$test.csv" \
    "$TX2GENE" "$ANALYSIS_ROOT/deseq2_$test" --test="$test" --replicates=3 "--ma-title=$title"
done
