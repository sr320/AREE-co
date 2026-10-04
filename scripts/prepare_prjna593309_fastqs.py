#!/usr/bin/env python3
"""Validate and stage the 12 PRJNA593309 paired-end libraries used in two AREE records.

Delisle et al. 2020 sequenced one whole oyster per tank (three tanks per temperature) at each
time after OsHV-1 cohabitation. AREE uses two 3 v 3 contrasts and downloads only their libraries:
OsHV-1 challenge at 21 C (48 h versus 0 h) and temperature during infection (29 C versus 21 C
at 12 h). The 26 C libraries and the other time points are not used.

``--download`` fetches each file with a single resumable connection. For a faster transfer,
fetch the manifest URLs with a multi-connection client into ``<analysis-root>/fastq`` as
``<run>_<mate>.fastq.gz.part`` and then run ``--download``: a complete ``.part`` file is
promoted only after it matches the manifest size and MD5. Use one connection per file.
"""

import argparse
import csv
import sys
from pathlib import Path

from aree.raw.fastq import download_runs, load_run_manifest, preflight, sample_name, write_design_samplesheet, write_nfcore_samplesheet


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OSHV1_TEMP_RNASEQ_PRJNA593309_runs.tsv")
EXPECTED_CONDITIONS = {"ctr_21C_0h": 3, "oshv_21C_48h": 3, "oshv_21C_12h": 3, "oshv_29C_12h": 3}
# Each record's design sheet relabels its two manifest groups as reference and test.
CONTRASTS = {
    "oshv1": {"control": "ctr_21C_0h", "oshv1": "oshv_21C_48h"},
    "heat": {"control": "oshv_21C_12h", "heat": "oshv_29C_12h"},
}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def contrast_rows(rows, contrast):
    """Design rows for one record; samples keep their manifest names so they match the Salmon output."""
    labels = {group: label for label, group in CONTRASTS[contrast].items()}
    return [{"sample": sample_name(row), "condition": labels[row["condition"]], "replicate": row["replicate"],
             "run_accession": row["run_accession"], "tank": "{}_tank{}".format(row["condition"].split("_")[1], row["replicate"])}
            for row in rows if row["condition"] in labels]


def write_contrast_sheet(rows, path):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "condition", "replicate", "run_accession", "tank"])
        writer.writeheader()
        writer.writerows(rows)
    print("wrote DESeq2 design sheet {}".format(path))


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
    print("runs: 12 paired-end (3 tanks x 21 C 0 h, 21 C 12 h, 21 C 48 h, 29 C 12 h)")
    if not preflight(rows, args.analysis_root, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    if args.download:
        download_runs(rows, output_dir, workers=args.workers, retries=args.retries, timeout=90)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # All 12 libraries, for the Salmon QC summary.
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv")
    for contrast in CONTRASTS:
        write_contrast_sheet(contrast_rows(rows, contrast), args.analysis_root / "deseq2_samplesheet_{}.csv".format(contrast))
    return 0


if __name__ == "__main__":
    sys.exit(main())
