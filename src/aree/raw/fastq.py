"""ENA run manifests, checksum-verified resumable FASTQ downloads, and sample sheets."""

import csv
import hashlib
import http.client
import os
import shutil
import socket
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


NETWORK_ERRORS = (socket.timeout, TimeoutError, ConnectionError, urllib.error.URLError, http.client.IncompleteRead)


def md5sum(path):
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_run_manifest(path, expected_conditions):
    """Read an ENA run manifest and check it has exactly the expected runs per condition.

    ``expected_conditions`` maps condition label to run count, e.g. ``{"control": 3, "heat": 3}``.
    """
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    counts = {}
    for row in rows:
        counts[row["condition"]] = counts.get(row["condition"], 0) + 1
    if counts != expected_conditions:
        raise ValueError("Expected runs per condition {}; found {}".format(expected_conditions, counts))
    if len({row["run_accession"] for row in rows}) != len(rows):
        raise ValueError("Run accessions must be unique")
    return rows


def total_fastq_bytes(rows):
    return sum(int(row[field]) for row in rows for field in ("fastq_1_bytes", "fastq_2_bytes"))


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


def preflight(rows, destination, scratch_multiplier):
    """Disk-space check for downloading a manifest's FASTQs into destination."""
    compressed = total_fastq_bytes(rows)
    required = int(compressed * scratch_multiplier)
    probe_dir = nearest_existing_ancestor(destination)
    free = shutil.disk_usage(str(probe_dir)).free
    print("compressed FASTQs: {:.2f} GB".format(compressed / 1e9))
    print("required free space: {:.2f} GB".format(required / 1e9))
    print("capacity checked at: {}".format(probe_dir))
    print("available free space: {:.2f} GB".format(free / 1e9))
    return free >= required


def download_once(url, destination, expected_bytes, expected_md5, timeout=60):
    """Download url to destination via a ``.part`` file, resuming it when the server allows.

    An existing destination or completed ``.part`` file is accepted only if its size and MD5
    match; nothing is ever promoted to destination without passing both checks.
    """
    destination = Path(destination)
    if destination.exists():
        if destination.stat().st_size == expected_bytes and md5sum(destination) == expected_md5:
            print("verified existing {}".format(destination), flush=True)
            return
        raise ValueError("Existing file fails size or MD5 validation: {}".format(destination))
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists() and partial.stat().st_size == expected_bytes:
        if md5sum(partial) == expected_md5:
            os.replace(str(partial), str(destination))
            print("verified completed partial {}".format(destination), flush=True)
            return
        raise ValueError("Completed partial file fails MD5 validation: {}".format(partial))
    if partial.exists() and partial.stat().st_size > expected_bytes:
        raise ValueError("Partial file is larger than expected: {}".format(partial))

    offset = partial.stat().st_size if partial.exists() else 0
    request = urllib.request.Request(url)
    if offset:
        request.add_header("Range", "bytes={}-".format(offset))
    with urllib.request.urlopen(request, timeout=timeout) as source:
        append = offset > 0 and getattr(source, "status", None) == 206
        if offset and not append:
            print("server did not honor resume request; restarting {}".format(destination), flush=True)
        elif append:
            print("resuming {} at byte {}".format(destination, offset), flush=True)
        with partial.open("ab" if append else "wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    if partial.stat().st_size != expected_bytes or md5sum(partial) != expected_md5:
        raise ValueError("Downloaded file fails size or MD5 validation: {}".format(partial))
    os.replace(str(partial), str(destination))
    print("downloaded and verified {}".format(destination), flush=True)


def download(url, destination, expected_bytes, expected_md5, retries=12, timeout=60):
    """``download_once`` with exponential backoff (capped at 30 s) on network errors."""
    for attempt in range(retries + 1):
        try:
            return download_once(url, destination, expected_bytes, expected_md5, timeout=timeout)
        except NETWORK_ERRORS as error:
            if attempt == retries:
                raise
            delay = min(2 ** attempt, 30)
            print(
                "network error for {} (attempt {}/{}): {}; retrying in {}s".format(
                    destination, attempt + 1, retries + 1, error, delay
                ),
                flush=True,
            )
            time.sleep(delay)


def fastq_path(output_dir, run_accession, mate):
    return Path(output_dir) / "{}_{}.fastq.gz".format(run_accession, mate)


def download_runs(rows, output_dir, workers=1, retries=12, timeout=60):
    """Download both mates of every manifest run into output_dir, verifying each file."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    jobs = [
        (
            row["fastq_{}".format(mate)],
            fastq_path(output_dir, row["run_accession"], mate),
            int(row["fastq_{}_bytes".format(mate)]),
            row["fastq_{}_md5".format(mate)],
        )
        for row in rows
        for mate in ("1", "2")
    ]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(download, *job, retries=retries, timeout=timeout) for job in jobs]
        for future in as_completed(futures):
            future.result()


def sample_name(row):
    return "{}_{}".format(row["condition"], row["replicate"])


def write_nfcore_samplesheet(rows, output_dir, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sample": sample_name(row),
                    "fastq_1": str(fastq_path(output_dir, row["run_accession"], "1").resolve()),
                    "fastq_2": str(fastq_path(output_dir, row["run_accession"], "2").resolve()),
                    "strandedness": "auto",
                }
            )
    print("wrote nf-core sample sheet {}".format(path))


def write_design_samplesheet(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "condition", "replicate", "run_accession"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sample": sample_name(row),
                    "condition": row["condition"],
                    "replicate": row["replicate"],
                    "run_accession": row["run_accession"],
                }
            )
    print("wrote DESeq2 design sheet {}".format(path))
