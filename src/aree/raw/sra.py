"""NCBI SRA archive verification and streaming FASTQ conversion helpers."""

import csv
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from aree.raw.fastq import fastq_path, md5sum


# Gzipped FASTQ bytes per spot, measured over both mates at fastp's compression level 4.
BYTES_PER_SPOT = 175


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


def check_free_space(rows, sra_dir, output_dir, keep_sra):
    """Return a complaint if the runs still to convert will not fit, else None."""
    pending = [
        row for row in rows
        if not all(fastq_path(output_dir, row["run_accession"], mate).exists() for mate in ("1", "2"))
    ]
    needed = sum(int(row["read_count"]) for row in pending) * BYTES_PER_SPOT
    if not keep_sra:
        # Each archive is removed once its FASTQ pair lands, so that space comes back.
        needed -= sum(
            path.stat().st_size
            for path in (sra_path(sra_dir, row["run_accession"]) for row in pending)
            if path.exists()
        )
    free = shutil.disk_usage(output_dir).free
    print("{} runs to convert, needing ~{:.0f} GB net; {:.0f} GB free".format(
        len(pending), needed / 1e9, free / 1e9))
    if needed > free * 0.9:
        return "not enough room: need ~{:.0f} GB, have {:.0f} GB free on {}".format(
            needed / 1e9, free / 1e9, output_dir)
    return None


def discard(paths):
    for path in paths:
        path.unlink(missing_ok=True)


def check_output(run, row, staged, report):
    """Reject output that exited cleanly but does not hold the run's full read complement.

    A full disk makes fastp lose writes without reporting failure, which leaves a file that
    still ends on the run's last spot and so survives any end-of-file check. Comparing the
    counts fastp itself reports against the manifest is what actually catches that.
    """
    for path in staged:
        if not path.exists():
            raise RuntimeError("conversion did not produce {}".format(path))
        with path.open("rb") as handle:
            if handle.read(2) != b"\x1f\x8b":
                raise RuntimeError("{} is not gzip; fastp wrote plain text".format(path))
    summary = json.loads(report.read_text())["summary"]["before_filtering"]
    expected = int(row["read_count"])
    # fastp counts both mates; the manifest counts spots.
    if summary["total_reads"] != expected * 2:
        raise RuntimeError(
            "{}: wrote {} reads, manifest expects {} ({} spots x 2)".format(
                run, summary["total_reads"], expected * 2, expected
            )
        )


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

    # Staged in a subdirectory so an interrupted run is not mistaken for a finished one.
    # The file names must be preserved: fastp chooses gzip output from the ``.gz`` suffix,
    # so renaming the staged copies to e.g. ``.fastq.gz.partial`` silently writes plain text.
    staging = Path(output_dir) / ".staging"
    staging.mkdir(exist_ok=True)
    staged = [staging / path.name for path in mates]
    report = staging / "{}.fastp.json".format(run)
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
            "--json", str(report),
            "--html", "/dev/null",
        ],
        stdin=dump.stdout,
    )
    dump.stdout.close()  # so fastq-dump sees EPIPE if fastp dies first
    split_code = split.wait()
    dump_code = dump.wait()
    if dump_code != 0 or split_code != 0:
        discard(staged + [report])
        raise RuntimeError(
            "{}: fastq-dump exited {}, fastp exited {}".format(run, dump_code, split_code)
        )
    try:
        check_output(run, row, staged, report)
    except Exception:
        discard(staged + [report])
        raise
    report.unlink(missing_ok=True)
    for source, target in zip(staged, mates):
        source.rename(target)
    if not keep_sra:
        archive.unlink()
    print("converted {}".format(run), flush=True)
