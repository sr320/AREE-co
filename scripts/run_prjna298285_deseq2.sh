#!/usr/bin/env bash
# DESeq2 for PRJNA298285: acidified (pH 7.9), warmed (22 C) and acidified + warmed larvae, each
# against the ambient larvae (pH 8.1, 20 C), blocking on developmental stage (trochophore, early and
# late veliger). Two libraries per stage and treatment.
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
for test in acidified warmed acidified_warmed; do
  case $test in
    acidified) title="pH 7.9 versus pH 8.1 at 20 C (larvae, stage-blocked)" ;;
    warmed) title="22 C versus 20 C at pH 8.1 (larvae, stage-blocked)" ;;
    acidified_warmed) title="pH 7.9 at 22 C versus pH 8.1 at 20 C (larvae, stage-blocked)" ;;
  esac
  Rscript "$SCRIPT_DIR/run_salmon_tximport_deseq2.R" "$QUANT_DIR" "$ANALYSIS_ROOT/deseq2_samplesheet_$test.csv" \
    "$TX2GENE" "$ANALYSIS_ROOT/deseq2_$test" --test="$test" --replicates=0 --covariate=stage "--ma-title=$title"
done
