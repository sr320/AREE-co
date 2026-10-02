#!/usr/bin/env python3
"""Validate and stage the 51 PRJNA735889 juvenile paired-end FASTQs used in the two-band pH contrast.

``--download`` fetches each file with a single resumable connection. For a faster transfer,
fetch the manifest URLs with a multi-connection client into ``<analysis-root>/fastq`` as
``<run>_<mate>.fastq.gz.part`` and then run ``--download``: a complete ``.part`` file is
promoted only after it matches the manifest size and MD5.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import download_runs, load_run_manifest, preflight, write_design_samplesheet, write_nfcore_samplesheet


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OA_RNASEQ_PRJNA735889_runs.tsv")
EXPECTED_CONDITIONS = {"control": 26, "acidified": 25}


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
    print("runs: 51 paired-end (26 at pHT 7.4-7.8, 25 at pHT 6.5-6.8; tipping-window tanks excluded)")
    if not preflight(rows, args.analysis_root, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    if args.download:
        download_runs(rows, output_dir, workers=args.workers, retries=args.retries, timeout=90)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # tank and pH are carried for QC; tank is nested in condition, so it is not a covariate.
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("tank", "ph_total"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
