"""Validate the PRJNA694496 manifest and optionally download its FASTQs."""

import argparse
import csv
import hashlib
import http.client
import os
import socket
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_MANIFEST = Path("data/manifests/CGIG_THERMOTOL_RNASEQ_PRJNA694496_runs.tsv")
DEFAULT_OUTPUT = Path("data/raw/CGIG_THERMOTOL_RNASEQ_PRJNA694496")


def load_manifest(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 6:
        raise ValueError("Expected six PRJNA694496 runs; found {}".format(len(rows)))
    if len({row["run_accession"] for row in rows}) != len(rows):
        raise ValueError("Run accessions must be unique")
    counts = {}
    for row in rows:
        counts[row["condition"]] = counts.get(row["condition"], 0) + 1
    if counts != {"selected": 3, "control": 3}:
        raise ValueError("Expected three selected and three control runs; found {}".format(counts))
    return rows


def total_fastq_bytes(rows):
    return sum(int(row[field]) for row in rows for field in ("fastq_1_bytes", "fastq_2_bytes"))


def md5sum(path):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def nearest_existing_ancestor(path):
    """Return the closest existing directory at or above path."""
    probe = Path(path).resolve()
    while not probe.exists():
        parent = probe.parent
        if parent == probe:
            raise FileNotFoundError("No existing ancestor for {}".format(path))
        probe = parent
    if not probe.is_dir():
        probe = probe.parent
    return probe


def download_once(url, destination, expected_bytes, expected_md5):
    destination = Path(destination)
    if destination.exists():
        if destination.stat().st_size == expected_bytes and md5sum(destination) == expected_md5:
            print("verified existing {}".format(destination))
            return
        raise ValueError("Existing file fails size or MD5 validation: {}".format(destination))
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists() and partial.stat().st_size == expected_bytes:
        if md5sum(partial) == expected_md5:
            os.replace(str(partial), str(destination))
            print("verified completed partial {}".format(destination))
            return
        raise ValueError("Completed partial file fails MD5 validation: {}".format(partial))
    if partial.exists() and partial.stat().st_size > expected_bytes:
        raise ValueError("Partial file is larger than expected: {}".format(partial))

    offset = partial.stat().st_size if partial.exists() else 0
    request = urllib.request.Request(url)
    if offset:
        request.add_header("Range", "bytes={}-".format(offset))
    with urllib.request.urlopen(request, timeout=60) as source:
        append = offset > 0 and getattr(source, "status", None) == 206
        mode = "ab" if append else "wb"
        if offset and not append:
            print("server did not honor resume request; restarting {}".format(destination))
        elif append:
            print("resuming {} at byte {}".format(destination, offset))
        with partial.open(mode) as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    if partial.stat().st_size != expected_bytes or md5sum(partial) != expected_md5:
        raise ValueError("Downloaded file fails size or MD5 validation: {}".format(partial))
    os.replace(str(partial), str(destination))
    print("downloaded and verified {}".format(destination))


def download(url, destination, expected_bytes, expected_md5, retries=12):
    network_errors = (
        socket.timeout,
        TimeoutError,
        ConnectionError,
        urllib.error.URLError,
        http.client.IncompleteRead,
    )
    for attempt in range(retries + 1):
        try:
            return download_once(url, destination, expected_bytes, expected_md5)
        except network_errors as error:
            if attempt == retries:
                raise
            delay = min(2 ** attempt, 30)
            print(
                "network error for {} (attempt {}/{}): {}; retrying in {}s".format(
                    destination, attempt + 1, retries + 1, error, delay
                )
            )
            time.sleep(delay)


def write_nfcore_samplesheet(rows, output_dir, path):
    fieldnames = ["sample", "fastq_1", "fastq_2", "strandedness"]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            run = row["run_accession"]
            writer.writerow(
                {
                    "sample": "{}_{}".format(row["condition"], row["replicate"]),
                    "fastq_1": str((output_dir / "{}_1.fastq.gz".format(run)).resolve()),
                    "fastq_2": str((output_dir / "{}_2.fastq.gz".format(run)).resolve()),
                    "strandedness": "auto",
                }
            )
    print("wrote nf-core sample sheet {}".format(path))


def write_design_samplesheet(rows, path):
    fieldnames = ["sample", "condition", "replicate", "run_accession"]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sample": "{}_{}".format(row["condition"], row["replicate"]),
                    "condition": row["condition"],
                    "replicate": row["replicate"],
                    "run_accession": row["run_accession"],
                }
            )
    print("wrote DESeq2 design sheet {}".format(path))


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
    compressed_bytes = total_fastq_bytes(rows)
    required_bytes = int(compressed_bytes * args.scratch_multiplier)
    probe_dir = nearest_existing_ancestor(args.output_dir)
    free_bytes = shutil.disk_usage(str(probe_dir)).free
    print("manifest runs: {} (3 selected, 3 control)".format(len(rows)))
    print("compressed FASTQs: {:.2f} GB".format(compressed_bytes / 1e9))
    print("required free space: {:.2f} GB".format(required_bytes / 1e9))
    print("capacity checked at: {}".format(probe_dir))
    print("available free space: {:.2f} GB".format(free_bytes / 1e9))
    if free_bytes < required_bytes:
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    print("preflight passed")
    if not args.download:
        if args.write_samplesheet:
            args.write_samplesheet.parent.mkdir(parents=True, exist_ok=True)
            write_nfcore_samplesheet(rows, args.output_dir, args.write_samplesheet)
        if args.write_design_sheet:
            args.write_design_sheet.parent.mkdir(parents=True, exist_ok=True)
            write_design_samplesheet(rows, args.write_design_sheet)
        return 0
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.retries < 0:
        parser.error("--retries cannot be negative")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for row in rows:
        run = row["run_accession"]
        for mate in ("1", "2"):
            jobs.append(
                (
                    row["fastq_{}".format(mate)],
                    args.output_dir / "{}_{}.fastq.gz".format(run, mate),
                    int(row["fastq_{}_bytes".format(mate)]),
                    row["fastq_{}_md5".format(mate)],
                    args.retries,
                )
            )
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(download, *job) for job in jobs]
        for future in as_completed(futures):
            future.result()
    if args.write_samplesheet:
        args.write_samplesheet.parent.mkdir(parents=True, exist_ok=True)
        write_nfcore_samplesheet(rows, args.output_dir, args.write_samplesheet)
    if args.write_design_sheet:
        args.write_design_sheet.parent.mkdir(parents=True, exist_ok=True)
        write_design_samplesheet(rows, args.write_design_sheet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
