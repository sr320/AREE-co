#!/usr/bin/env python3
"""Stage the 42 PRJNA1329250 control/OsHV-1 FASTQ pairs from NCBI SRA archives.

Unlike the ENA-sourced studies, this BioProject is served far faster from the NCBI
Open Data S3 mirror, which publishes normalized ``.sra`` archives rather than the
submitter's FASTQ files. Acquisition is therefore two stages: verify the ``.sra``
archives against the MD5s NCBI publishes for them, then convert each archive to a
FASTQ pair. The run manifest still records the ENA FASTQ URLs and checksums, but those
describe the mirror we did not use, so this script validates against
``*_sra_locations.tsv`` instead.

Conversion streams rather than using ``fasterq-dump``. ``fasterq-dump`` writes roughly
20x the archive size to a scratch directory before emitting anything, which on the
external volume this study is staged on costs about 45 minutes per run. Legacy
``fastq-dump --stdout`` uses no scratch at all and emits the mates interleaved, so the
only remaining work is splitting and compressing that stream; ``fastp`` does both with
all of its filtering disabled, at close to the speed of the extraction itself.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import load_run_manifest, write_design_samplesheet, write_nfcore_samplesheet
from aree.raw.sra import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OSHV1_RNASEQ_PRJNA1329250_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_OSHV1_RNASEQ_PRJNA1329250_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 12, "oshv1_uvar": 30}


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
    # The archives live on a single external disk, so concurrent readers only cause seeking.
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
    print("runs: 42 (12 control, 30 oshv1_uvar)")

    missing_locations = sorted({row["run_accession"] for row in rows} - set(locations))
    if missing_locations:
        print("no SRA location recorded for: {}".format(missing_locations), file=sys.stderr)
        return 2

    if args.verify:
        failures = verify_sra(rows, locations, sra_dir, args.workers)
        if failures:
            print("SRA verification failed for: {}".format(sorted(failures)), file=sys.stderr)
            return 2
        print("all 42 SRA archives match the NCBI size and MD5")

    if args.convert:
        shortfall = check_free_space(rows, sra_dir, output_dir, args.keep_sra)
        if shortfall:
            print(shortfall, file=sys.stderr)
            return 2
        # Conversion is disk-bound; running it serially keeps the external volume sequential.
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            convert_run(
                row, sra_dir, output_dir, args.threads,
                args.fastq_dump, args.fastp, args.keep_sra,
            )

    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # population is carried through so DESeq2 can block on oyster lineage.
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("population",))
    return 0


if __name__ == "__main__":
    sys.exit(main())
