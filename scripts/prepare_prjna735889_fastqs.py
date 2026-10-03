#!/usr/bin/env python3
"""Validate and stage the 51 PRJNA735889 juvenile paired-end runs used in the two-band pH contrast.

All runs come from the NCBI Open Data S3 mirror as ``.sra`` archives listed in
``*_sra_locations.tsv``: verified against NCBI's published size and MD5 (``--verify``), then
streamed to gzipped FASTQ pairs whose read counts must match the manifest (``--convert``),
using the steps shared with ``prepare_prjna1329250_fastqs``. ENA served this study at about
80 KB/s on 2 October 2026 (about a month for 238 GB), so it was not used.

Because it imports the shared SRA steps from ``scripts``, run it from the repository root
as ``python -m scripts.prepare_prjna735889_fastqs``.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import load_run_manifest, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OA_RNASEQ_PRJNA735889_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_OA_RNASEQ_PRJNA735889_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 26, "acidified": 25}
# 13Ebis is a second library from oyster 13E (expressed-SNP allele frequencies differ by 0.026,
# below the 0.032 expected from resequencing one RNA), so it is kept in the manifest for
# provenance but left out of the sample sheets to avoid counting one animal twice.
EXCLUDED_RUNS = {"SRR14803132"}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


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
    print("runs: 51 paired-end (26 at pHT 7.4-7.8, 25 at pHT 6.5-6.8; tipping-window tanks excluded); 50 analysed")
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
    analysed = [row for row in rows if row["run_accession"] not in EXCLUDED_RUNS]
    write_nfcore_samplesheet(analysed, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # tank and pH are carried for QC; tank is nested in condition, so it is not a covariate.
    write_design_samplesheet(analysed, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("tank", "ph_total"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
