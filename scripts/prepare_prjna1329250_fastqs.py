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


def convert_run(row, sra_dir, output_dir, threads, fastq_dump, fastp, keep_sra):
    """Convert one ``.sra`` archive to a gzipped FASTQ pair, skipping work already done."""
    run = row["run_accession"]
    mates = [fastq_path(output_dir, run, mate) for mate in ("1", "2")]
    if all(path.exists() for path in mates):
        print("skip {} (already converted)".format(run), flush=True)
        return
    archive = sra_path(sra_dir, run)
    if not archive.exists():
        raise FileNotFoundError("Missing SRA archive for {}: {}".format(run, archive))

    # Written under temporary names so an interrupted run is not mistaken for a finished one.
    staged = [path.with_suffix(path.suffix + ".partial") for path in mates]
    dump = subprocess.Popen(
        [fastq_dump, "--split-spot", "--skip-technical", "--stdout", str(archive)],
        stdout=subprocess.PIPE,
    )
    split = subprocess.Popen(
        [
            fastp,
            "--stdin", "--interleaved_in",
            "--out1", str(staged[0]),
            "--out2", str(staged[1]),
            # Salmon is given untrimmed reads, so fastp is used purely as a splitter.
            "--disable_adapter_trimming",
            "--disable_quality_filtering",
            "--disable_length_filtering",
            "--disable_trim_poly_g",
            "--compression", "4",
            "--thread", str(threads),
            "--json", "/dev/null",
            "--html", "/dev/null",
        ],
        stdin=dump.stdout,
    )
    dump.stdout.close()  # so fastq-dump sees EPIPE if fastp dies first
    split_code = split.wait()
    dump_code = dump.wait()
    if dump_code != 0 or split_code != 0:
        for path in staged:
            path.unlink(missing_ok=True)
        raise RuntimeError(
            "{}: fastq-dump exited {}, fastp exited {}".format(run, dump_code, split_code)
        )
    for source, target in zip(staged, mates):
        if not source.exists():
            raise RuntimeError("conversion did not produce {}".format(target))
        source.rename(target)
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
