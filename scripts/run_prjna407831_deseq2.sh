#!/usr/bin/env bash
# DESeq2 for PRJNA407831: 35 C seawater for 6 h and for 24 h, each against 0 h, blocking on the four
# source populations. Each library pools gill RNA of five oysters from one tank, and the same tanks
# are sampled at every time, so a tank-paired model is kept as a sensitivity analysis.
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
for test in h6 h24; do
  hours=${test#h}
  Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_$test.csv" \
    "$TX2GENE" "$ANALYSIS_ROOT/deseq2_$test" --test="$test" --replicates=0 --covariate=population \
    "--ma-title=35 C seawater for $hours h versus 0 h (gill)"
  Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_$test.csv" \
    "$TX2GENE" "$ANALYSIS_ROOT/deseq2_${test}_tank_paired" --test="$test" --replicates=0 --covariate=tank \
    "--ma-title=35 C seawater for $hours h versus 0 h (gill), tank-paired"
done
