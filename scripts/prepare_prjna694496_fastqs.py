"""Validate the PRJNA694496 manifest and optionally download its FASTQs."""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import (  # noqa: F401  (re-exported for callers of this script's module)
    download,
    download_runs,
    load_run_manifest,
    md5sum,
    nearest_existing_ancestor,
    preflight,
    total_fastq_bytes,
    write_design_samplesheet,
    write_nfcore_samplesheet,
)


DEFAULT_MANIFEST = Path("data/manifests/CGIG_THERMOTOL_RNASEQ_PRJNA694496_runs.tsv")
DEFAULT_OUTPUT = Path("data/raw/CGIG_THERMOTOL_RNASEQ_PRJNA694496")
EXPECTED_CONDITIONS = {"selected": 3, "control": 3}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def write_sheets(args, rows):
    if args.write_samplesheet:
        write_nfcore_samplesheet(rows, args.output_dir, args.write_samplesheet)
    if args.write_design_sheet:
        write_design_samplesheet(rows, args.write_design_sheet)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--scratch-multiplier", type=float, default=1.85)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--retries", type=int, default=12)
    parser.add_argument("--write-samplesheet", type=Path)
    parser.add_argument("--write-design-sheet", type=Path)
    args = parser.parse_args()

    rows = load_manifest(args.manifest)
    print("manifest runs: {} (3 selected, 3 control)".format(len(rows)))
    if not preflight(rows, args.output_dir, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    print("preflight passed")
    if args.download:
        if args.workers < 1:
            parser.error("--workers must be at least 1")
        if args.retries < 0:
            parser.error("--retries cannot be negative")
        download_runs(rows, args.output_dir, workers=args.workers, retries=args.retries, timeout=60)
    write_sheets(args, rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
