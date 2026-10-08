#!/usr/bin/env bash
# DESeq2 for PRJEB58194 (PESTO): three records, each blocking on developmental stage.
#   f0_exposed:  F0 embryos and larvae exposed to the pesticide mixture against F0 controls
#   f1_direct:   F1 offspring exposed (TE) against F1 controls (TT)
#   f1_parental: F1 control offspring of exposed parents (ET) against TT
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
for test in f0_exposed f1_direct f1_parental; do
  case $test in
    f0_exposed) title="F0 pesticide mixture 0-48 hpf versus control (gastrula and D larva)" ;;
    f1_direct) title="F1 exposed offspring of control parents (TE) versus TT" ;;
    f1_parental) title="F1 control offspring of exposed parents (ET) versus TT" ;;
  esac
  Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_$test.csv" \
    "$TX2GENE" "$ANALYSIS_ROOT/deseq2_$test" --test="$test" --replicates=0 --covariate=stage "--ma-title=$title"
done
