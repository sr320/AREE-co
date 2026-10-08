#!/usr/bin/env python3
"""Validate and stage the 24 PRJNA298285 larval paired-end libraries (pH x temperature x stage).

Larvae from Point du Chateau (Brittany) were reared at pH 8.1 or 7.9 and 20 or 22 C and sampled as
trochophores (14 h), early veligers (24 h) and late veligers (36 h), two libraries per cell. AREE
tests acidified, warmed and acidified + warmed larvae each against the ambient larvae, blocking
on stage.

All runs come from the NCBI Open Data S3 mirror as ``.sra`` archives listed in
``*_sra_locations.tsv``: verified against NCBI's published size and MD5 (``--verify``), then
streamed to gzipped FASTQ pairs whose read counts must match the manifest (``--convert``),
using the steps shared with ``prepare_prjna1329250_fastqs``.

Because it imports the shared SRA steps from ``scripts``, run it from the repository root
as ``python -m scripts.prepare_prjna298285_fastqs``.
"""

import argparse
import csv
import sys
from pathlib import Path

from aree.raw.fastq import load_run_manifest, sample_name, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_LARVAL_OA_TEMP_RNASEQ_PRJNA298285_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_LARVAL_OA_TEMP_RNASEQ_PRJNA298285_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 6, "acidified": 6, "warmed": 6, "acidified_warmed": 6}
# Each treatment is tested against the ambient larvae in its own design sheet.
CONTRASTS = ("acidified", "warmed", "acidified_warmed")
# Libraries left out of the DESeq2 design sheets after QC (none yet).
EXCLUDED_SAMPLES = set()


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def write_contrast_sheet(rows, test, path):
    """One treatment against ambient larvae; samples keep their manifest names for the Salmon output."""
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "condition", "replicate", "run_accession", "stage"])
        writer.writeheader()
        for row in rows:
            if row["condition"] in ("control", test) and sample_name(row) not in EXCLUDED_SAMPLES:
                writer.writerow({"sample": sample_name(row), "condition": row["condition"],
                                 "replicate": row["replicate"], "run_accession": row["run_accession"], "stage": row["stage"]})
    print("wrote DESeq2 design sheet {}".format(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locations", type=Path, default=DEFAULT_LOCATIONS)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--sra-dir", type=Path, help="defaults to <analysis-root>/sra")
    parser.add_argument("--verify", action="store_true", help="check the .sra archives against NCBI MD5s")
    parser.add_argument("--convert", action="store_true", help="stream each .sra to a gzipped FASTQ pair")
    parser.add_argument("--keep-sra", action="store_true", help="retain .sra archives after conversion")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--threads", type=int, default=8, help="fastp compression threads")
    parser.add_argument("--fastq-dump", default="fastq-dump")
    parser.add_argument("--fastp", default="fastp")
    args = parser.parse_args()
    rows = load_manifest(args.manifest)
    locations = load_locations(args.locations)
    if set(locations) != {row["run_accession"] for row in rows}:
        raise ValueError("SRA locations and manifest name different runs")
    output_dir = args.analysis_root / "fastq"
    output_dir.mkdir(parents=True, exist_ok=True)
    sra_dir = args.sra_dir or args.analysis_root / "sra"
    print("runs: 24 paired-end larval libraries (2 pH x 2 temperatures x 3 stages x 2)")
    if args.verify:
        failures = verify_sra(rows, locations, sra_dir, args.workers)
        if failures:
            print("SRA verification failed for: {}".format(sorted(failures)), file=sys.stderr)
            return 2
        print("all {} SRA archives match the NCBI size and MD5".format(len(rows)))
    if args.convert:
        shortfall = check_free_space(rows, sra_dir, output_dir, args.keep_sra)
        if shortfall:
            print(shortfall, file=sys.stderr)
            return 2
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            convert_run(row, sra_dir, output_dir, args.threads, args.fastq_dump, args.fastp, args.keep_sra)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("stage",))
    for test in CONTRASTS:
        write_contrast_sheet(rows, test, args.analysis_root / "deseq2_samplesheet_{}.csv".format(test))
    return 0


if __name__ == "__main__":
    sys.exit(main())
