#!/usr/bin/env python3
"""Validate and stage the 36 PRJNA407831 gill paired-end libraries (Li et al. 2018 heat stress).

Each library pools gill RNA of five F1 oysters from one culture tank (four source populations x
three tanks), sampled after 0, 6 and 24 h in 35 C seawater; the same tanks are sampled at each
time. AREE tests 6 h and 24 h each against 0 h, blocking on population.

All runs come from the NCBI Open Data S3 mirror as ``.sra`` archives listed in
``*_sra_locations.tsv``: verified against NCBI's published size and MD5 (``--verify``), then
streamed to gzipped FASTQ pairs whose read counts must match the manifest (``--convert``),
using the steps shared with ``prepare_prjna1329250_fastqs``. ENA served other studies at
0.05-0.35 MB/s in October 2026 (1-3 days for these 94 GB), so it was not used.

Because it imports the shared SRA steps from ``scripts``, run it from the repository root
as ``python -m scripts.prepare_prjna407831_fastqs``.
"""

import argparse
import csv
import sys
from pathlib import Path

from aree.raw.fastq import load_run_manifest, sample_name, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_HEAT_RNASEQ_PRJNA407831_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_HEAT_RNASEQ_PRJNA407831_sra_locations.tsv")
EXPECTED_CONDITIONS = {"h0": 12, "h6": 12, "h24": 12}
# Each heat time is tested against 0 h in its own design sheet.
CONTRASTS = ("h6", "h24")
def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def write_contrast_sheet(rows, test, path):
    """One heat time against 0 h; samples keep their manifest names so they match the Salmon output."""
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "condition", "replicate", "run_accession", "population", "tank"])
        writer.writeheader()
        for row in rows:
            if row["condition"] in ("h0", test):
                writer.writerow({"sample": sample_name(row), "condition": "control" if row["condition"] == "h0" else test,
                                 "replicate": row["replicate"], "run_accession": row["run_accession"],
                                 "population": row["population"], "tank": row["tank"]})
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
    print("runs: 36 paired-end gill pools (4 populations x 3 tanks x 0, 6, 24 h at 35 C)")
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
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("population", "tank"))
    for test in CONTRASTS:
        write_contrast_sheet(rows, test, args.analysis_root / "deseq2_samplesheet_{}.csv".format(test))
    return 0


if __name__ == "__main__":
    sys.exit(main())
