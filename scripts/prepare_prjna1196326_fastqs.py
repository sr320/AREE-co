#!/usr/bin/env python3
"""Validate and stage the 6 PRJNA1196326 eyed-larva paired-end FASTQs from ENA.

``--download`` fetches each file with a single resumable connection. For a faster transfer,
fetch the manifest URLs with a multi-connection client into ``<analysis-root>/fastq`` as
``<run>_<mate>.fastq.gz.part`` and then run ``--download``: a complete ``.part`` file is
promoted only after it matches the manifest size and MD5.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import download_runs, load_run_manifest, preflight, write_design_samplesheet, write_nfcore_samplesheet


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OA_RNASEQ_PRJNA1196326_runs.tsv")
EXPECTED_CONDITIONS = {"control": 3, "parental_acidified": 3}


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
    print("runs: 6 paired-end (3 ambient-parent S81, 3 acidified-parent S74)")
    if not preflight(rows, args.analysis_root, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    if args.download:
        download_runs(rows, output_dir, workers=args.workers, retries=args.retries, timeout=90)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
