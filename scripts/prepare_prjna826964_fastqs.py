#!/usr/bin/env python3
"""Stage the 18 PRJNA826964 control/acidified hepatopancreas FASTQ pairs from NCBI SRA archives.

The ENA FASTQ mirror served this BioProject at well under 100 KB/s, so, as for
PRJNA1329250, reads are taken from the NCBI Open Data S3 mirror as ``.sra`` archives,
verified against the size and MD5 NCBI publishes (``*_sra_locations.tsv``), and streamed
to gzipped FASTQ pairs with ``fastq-dump --stdout | fastp``. The verification and
conversion steps are shared with ``prepare_prjna1329250_fastqs``; the run manifest keeps
the ENA URLs and checksums for reference only.

Because it imports those shared steps from ``scripts``, run it from the repository root
as ``python -m scripts.prepare_prjna826964_fastqs``.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import load_run_manifest, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OA_RNASEQ_PRJNA826964_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_OA_RNASEQ_PRJNA826964_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 9, "acidified": 9}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locations", type=Path, default=DEFAULT_LOCATIONS)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--sra-dir", type=Path, help="defaults to <analysis-root>/sra")
    parser.add_argument("--verify", action="store_true", help="check .sra archives against NCBI MD5s")
    parser.add_argument("--convert", action="store_true", help="stream each .sra to a gzipped FASTQ pair")
    parser.add_argument("--keep-sra", action="store_true", help="retain .sra archives after conversion")
    parser.add_argument("--workers", type=int, default=1, help="parallel MD5 verifications")
    parser.add_argument("--threads", type=int, default=8, help="fastp compression threads")
    parser.add_argument("--fastq-dump", default="fastq-dump")
    parser.add_argument("--fastp", default="fastp")
    args = parser.parse_args()

    rows = load_manifest(args.manifest)
    locations = load_locations(args.locations)
    sra_dir = args.sra_dir or args.analysis_root / "sra"
    output_dir = args.analysis_root / "fastq"
    output_dir.mkdir(parents=True, exist_ok=True)
    print("runs: 18 (9 control, 9 acidified; 3 per arm at 7, 28 and 56 days)")

    missing_locations = sorted({row["run_accession"] for row in rows} - set(locations))
    if missing_locations:
        print("no SRA location recorded for: {}".format(missing_locations), file=sys.stderr)
        return 2

    if args.verify:
        failures = verify_sra(rows, locations, sra_dir, args.workers)
        if failures:
            print("SRA verification failed for: {}".format(sorted(failures)), file=sys.stderr)
            return 2
        print("all 18 SRA archives match the NCBI size and MD5")

    if args.convert:
        shortfall = check_free_space(rows, sra_dir, output_dir, args.keep_sra)
        if shortfall:
            print(shortfall, file=sys.stderr)
            return 2
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            convert_run(
                row, sra_dir, output_dir, args.threads,
                args.fastq_dump, args.fastp, args.keep_sra,
            )

    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # timepoint is carried so DESeq2 can block on exposure duration.
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("timepoint",))
    return 0


if __name__ == "__main__":
    sys.exit(main())
