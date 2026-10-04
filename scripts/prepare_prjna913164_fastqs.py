#!/usr/bin/env python3
"""Validate and stage the 72 PRJNA913164 QuantSeq 3' tag-seq single-end FASTQs from ENA.

``--download`` fetches each file with a single resumable connection. For a faster transfer,
fetch the manifest URLs with a multi-connection client into ``<analysis-root>/fastq`` as
``<run>_1.fastq.gz.part`` and then run ``--download``: a complete ``.part`` file is
promoted only after it matches the manifest size and MD5. Use one connection per file.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import download_runs, load_run_manifest, preflight, write_design_samplesheet, write_nfcore_samplesheet


DEFAULT_MANIFEST = Path("data/manifests/CGIG_HEATWAVE_RNASEQ_PRJNA913164_runs.tsv")
EXPECTED_CONDITIONS = {"control": 24, "heat": 24, "heat_emersion": 24}
# Each stressor is tested against the shared 20 C controls in its own design sheet.
CONTRASTS = ("heat", "heat_emersion")
# Libraries left out of the DESeq2 design sheets after QC: median VST correlation with the
# other libraries below the median minus three robust SDs (0.834) in a ~ ploidy + condition run on
# all 72. Seven match the authors' exclusions; their eighth, D54, passes and is kept. R62 and N57
# fall below the same cut-off and are also left out. See the study YAMLs.
EXCLUDED_SAMPLES = {"N56", "X44", "X42", "T62", "R53", "N54", "M43", "R62", "N57"}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--retries", type=int, default=12)
    parser.add_argument("--scratch-multiplier", type=float, default=1.25)
    args = parser.parse_args()
    rows = load_manifest(args.manifest)
    output_dir = args.analysis_root / "fastq"
    print("runs: 72 single-end tag-seq (control, 30 C, 30 C + 44 C emersion; 12 diploid and 12 triploid each)")
    if not preflight(rows, args.analysis_root, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    if args.download:
        download_runs(rows, output_dir, workers=args.workers, retries=args.retries, timeout=90)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # ploidy is carried so DESeq2 can block on it.
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("ploidy",))
    analysed = [row for row in rows if row["replicate"] not in EXCLUDED_SAMPLES]
    for test in CONTRASTS:
        subset = [row for row in analysed if row["condition"] in ("control", test)]
        write_design_samplesheet(subset, args.analysis_root / "deseq2_samplesheet_{}.csv".format(test), extra_columns=("ploidy",))
    return 0


if __name__ == "__main__":
    sys.exit(main())
