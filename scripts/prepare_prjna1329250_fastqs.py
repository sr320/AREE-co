#!/usr/bin/env python3
"""Stage the 42 PRJNA1329250 control/OsHV-1 FASTQ pairs from NCBI SRA archives.

Unlike the ENA-sourced studies, this BioProject is served far faster from the NCBI
Open Data S3 mirror, which publishes normalized ``.sra`` archives rather than the
submitter's FASTQ files. Acquisition is therefore two stages: verify the ``.sra``
archives against the MD5s NCBI publishes for them, then convert each archive to a
FASTQ pair with ``fasterq-dump``. The run manifest still records the ENA FASTQ URLs
and checksums, but those describe the mirror we did not use, so this script validates
against ``*_sra_locations.tsv`` instead.
"""

import argparse
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from aree.raw.fastq import fastq_path, load_run_manifest, md5sum, write_design_samplesheet, write_nfcore_samplesheet


DEFAULT_MANIFEST = Path("data/manifests/CGIG_OSHV1_RNASEQ_PRJNA1329250_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_OSHV1_RNASEQ_PRJNA1329250_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 12, "oshv1_uvar": 30}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def load_locations(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return {row["run_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}


def sra_path(sra_dir, run_accession):
    return Path(sra_dir) / "{}.sra".format(run_accession)


def verify_sra(rows, locations, sra_dir, workers):
    """Check every ``.sra`` archive against the size and MD5 that NCBI publishes for it."""

    def check(row):
        run = row["run_accession"]
        expected = locations[run]
        path = sra_path(sra_dir, run)
        if not path.exists():
            return run, "MISSING"
        size = path.stat().st_size
        if size != int(expected["sra_size"]):
            return run, "SIZE {} != {}".format(size, expected["sra_size"])
        if md5sum(path) != expected["sra_md5"]:
            return run, "MD5 MISMATCH"
        return run, "OK"

    failures = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(check, row): row for row in rows}
        for future in as_completed(futures):
            run, status = future.result()
            print("{} {}".format(status, run), flush=True)
            if status != "OK":
                failures.append(run)
    return failures


def convert_run(row, sra_dir, output_dir, temp_dir, threads, fasterq_dump, pigz, keep_sra):
    """Convert one ``.sra`` archive to a gzipped FASTQ pair, skipping work already done."""
    run = row["run_accession"]
    mates = [fastq_path(output_dir, run, mate) for mate in ("1", "2")]
    if all(path.exists() for path in mates):
        print("skip {} (already converted)".format(run), flush=True)
        return
    archive = sra_path(sra_dir, run)
    if not archive.exists():
        raise FileNotFoundError("Missing SRA archive for {}: {}".format(run, archive))

    subprocess.run(
        [
            fasterq_dump,
            "--split-files",
            "--skip-technical",
            "--threads", str(threads),
            "--temp", str(temp_dir),
            "--outdir", str(output_dir),
            str(archive),
        ],
        check=True,
    )
    plain = [Path(output_dir) / "{}_{}.fastq".format(run, mate) for mate in ("1", "2")]
    missing = [path for path in plain if not path.exists()]
    if missing:
        raise RuntimeError("fasterq-dump did not produce {}".format(missing))
    subprocess.run([pigz, "-p", str(threads), "--force", *[str(path) for path in plain]], check=True)
    for path in mates:
        if not path.exists():
            raise RuntimeError("compression did not produce {}".format(path))
    if not keep_sra:
        archive.unlink()
    print("converted {}".format(run), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locations", type=Path, default=DEFAULT_LOCATIONS)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--sra-dir", type=Path, help="defaults to <analysis-root>/sra")
    parser.add_argument("--verify", action="store_true", help="check .sra archives against NCBI MD5s")
    parser.add_argument("--convert", action="store_true", help="run fasterq-dump and compress the output")
    parser.add_argument("--keep-sra", action="store_true", help="retain .sra archives after conversion")
    parser.add_argument("--workers", type=int, default=4, help="parallel MD5 verifications")
    parser.add_argument("--threads", type=int, default=8, help="threads per fasterq-dump/pigz call")
    parser.add_argument("--fasterq-dump", default="fasterq-dump")
    parser.add_argument("--pigz", default="pigz")
    args = parser.parse_args()

    rows = load_manifest(args.manifest)
    locations = load_locations(args.locations)
    sra_dir = args.sra_dir or args.analysis_root / "sra"
    output_dir = args.analysis_root / "fastq"
    temp_dir = args.analysis_root / "tmp"
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
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
        # Conversion is disk-bound; running it serially keeps peak scratch use to one run.
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            convert_run(
                row, sra_dir, output_dir, temp_dir, args.threads,
                args.fasterq_dump, args.pigz, args.keep_sra,
            )

    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # population is carried through so DESeq2 can block on oyster lineage.
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("population",))
    return 0


if __name__ == "__main__":
    sys.exit(main())
