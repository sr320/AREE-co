#!/usr/bin/env python3
"""Validate and download the six PRJNA516762 control/heat FASTQ pairs."""

import argparse
import csv
import hashlib
import http.client
import os
import shutil
import socket
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_MANIFEST = Path("data/manifests/CGIG_HEAT_RNASEQ_PRJNA516762_runs.tsv")


def md5sum(path):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    counts = {}
    for row in rows:
        counts[row["condition"]] = counts.get(row["condition"], 0) + 1
    if len(rows) != 6 or counts != {"control": 3, "heat": 3}:
        raise ValueError("Expected three control and three heat runs")
    if len({row["run_accession"] for row in rows}) != 6:
        raise ValueError("Run accessions must be unique")
    return rows


def nearest_existing_ancestor(path):
    probe = Path(path).resolve()
    while not probe.exists():
        probe = probe.parent
    return probe if probe.is_dir() else probe.parent


def download_once(url, destination, expected_bytes, expected_md5):
    destination = Path(destination)
    partial = destination.with_suffix(destination.suffix + ".part")
    if destination.exists():
        if destination.stat().st_size == expected_bytes and md5sum(destination) == expected_md5:
            print("verified existing {}".format(destination), flush=True)
            return
        raise ValueError("Existing file fails size or MD5 validation: {}".format(destination))
    if partial.exists() and partial.stat().st_size == expected_bytes:
        if md5sum(partial) == expected_md5:
            os.replace(str(partial), str(destination))
            print("verified completed partial {}".format(destination), flush=True)
            return
        raise ValueError("Completed partial fails MD5 validation: {}".format(partial))
    offset = partial.stat().st_size if partial.exists() else 0
    if offset > expected_bytes:
        raise ValueError("Partial file is larger than expected: {}".format(partial))
    request = urllib.request.Request(url)
    if offset:
        request.add_header("Range", "bytes={}-".format(offset))
    with urllib.request.urlopen(request, timeout=90) as source:
        append = offset > 0 and getattr(source, "status", None) == 206
        if offset and not append:
            print("server did not honor resume; restarting {}".format(destination), flush=True)
        with partial.open("ab" if append else "wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    if partial.stat().st_size != expected_bytes or md5sum(partial) != expected_md5:
        raise ValueError("Downloaded file fails size or MD5 validation: {}".format(partial))
    os.replace(str(partial), str(destination))
    print("downloaded and verified {}".format(destination), flush=True)


def download(url, destination, expected_bytes, expected_md5, retries):
    errors = (socket.timeout, TimeoutError, ConnectionError, urllib.error.URLError,
              http.client.IncompleteRead)
    for attempt in range(retries + 1):
        try:
            return download_once(url, destination, expected_bytes, expected_md5)
        except errors as error:
            if attempt == retries:
                raise
            delay = min(2 ** attempt, 30)
            print("network error for {}: {}; retrying in {}s".format(
                destination, error, delay), flush=True)
            time.sleep(delay)


def write_sheets(rows, output_dir, nfcore_path, design_path):
    nfcore_path.parent.mkdir(parents=True, exist_ok=True)
    with nfcore_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        for row in rows:
            sample = "{}_{}".format(row["condition"], row["replicate"])
            writer.writerow({
                "sample": sample,
                "fastq_1": str((output_dir / (row["run_accession"] + "_1.fastq.gz")).resolve()),
                "fastq_2": str((output_dir / (row["run_accession"] + "_2.fastq.gz")).resolve()),
                "strandedness": "auto",
            })
    with design_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "condition", "replicate", "run_accession"])
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "sample": "{}_{}".format(row["condition"], row["replicate"]),
                "condition": row["condition"],
                "replicate": row["replicate"],
                "run_accession": row["run_accession"],
            })


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
    compressed = sum(int(row[key]) for row in rows for key in ("fastq_1_bytes", "fastq_2_bytes"))
    free = shutil.disk_usage(str(nearest_existing_ancestor(args.analysis_root))).free
    required = int(compressed * args.scratch_multiplier)
    print("runs: 6 (3 control, 3 heat)")
    print("compressed FASTQs: {:.2f} GB".format(compressed / 1e9))
    print("required free space: {:.2f} GB; available: {:.2f} GB".format(required / 1e9, free / 1e9))
    if free < required:
        return 2
    if args.download:
        output_dir.mkdir(parents=True, exist_ok=True)
        jobs = []
        for row in rows:
            for mate in ("1", "2"):
                jobs.append((row["fastq_" + mate],
                             output_dir / "{}_{}.fastq.gz".format(row["run_accession"], mate),
                             int(row["fastq_{}_bytes".format(mate)]),
                             row["fastq_{}_md5".format(mate)], args.retries))
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(download, *job) for job in jobs]
            for future in as_completed(futures):
                future.result()
    write_sheets(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv",
                 args.analysis_root / "deseq2_samplesheet.csv")
    print("wrote sample sheets under {}".format(args.analysis_root))
    return 0


if __name__ == "__main__":
    sys.exit(main())
