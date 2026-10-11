#!/usr/bin/env bash
# DESeq2 for PRJNA1244231 (gill): seawater-injected oysters at 12 h against 0 h, 3 v 3. "0" and "12"
# are read as hours of heat stress. The three libraries per group are one RNA pool resequenced, so
# this run documents the not-harmonized decision rather than producing a record.
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
Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_sw_12h.csv" \
  "$TX2GENE" "$ANALYSIS_ROOT/deseq2_sw_12h" --test=sw_12h --replicates=3 \
  "--ma-title=Seawater-injected 12 h versus 0 h (gill)"
