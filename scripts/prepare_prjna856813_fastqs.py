#!/usr/bin/env python3
"""Validate and stage the 18 PRJNA856813 gill single-end libraries from NCBI .sra archives.

AREE uses two day-41 gill contrasts from Bernardini et al. 2024: high-concentration tire leachate
and high-concentration tire microparticles, each against the shared controls (6 oysters per arm).

The ENA portal was unavailable when this study was staged (October 2026), so its FASTQ MD5s could
not be read. The runs come instead from the NCBI Open Data S3 mirror as ``.sra`` archives listed in
``*_sra_locations.tsv``: verified against NCBI's published size and MD5 (``--verify``), then streamed
to gzipped single-end FASTQ whose read count must match the SRA spot count (``--convert``).

Because it imports the shared SRA steps from ``scripts``, run it from the repository root as
``python -m scripts.prepare_prjna856813_fastqs``.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from aree.raw.fastq import fastq_path, load_run_manifest, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import discard, load_locations, sra_path, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_TIRE_RNASEQ_PRJNA856813_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_TIRE_RNASEQ_PRJNA856813_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 6, "leachate_high": 6, "tire_particles_high": 6}
# Each exposure is tested against the shared gill controls in its own design sheet.
CONTRASTS = ("leachate_high", "tire_particles_high")


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def convert_single_end(row, sra_dir, output_dir, threads, fastq_dump, fastp):
    """Stream one single-end ``.sra`` archive to gzipped FASTQ, checking the read count."""
    run = row["run_accession"]
    target = fastq_path(output_dir, run, "1")
    if target.exists():
        print("skip {} (already converted)".format(run), flush=True)
        return
    archive = sra_path(sra_dir, run)
    if not archive.exists():
        raise FileNotFoundError("Missing SRA archive for {}: {}".format(run, archive))
    # Staged under the real name in a subdirectory: fastp picks gzip output from the .gz suffix.
    staging = Path(output_dir) / ".staging"
    staging.mkdir(exist_ok=True)
    staged = staging / target.name
    report = staging / "{}.fastp.json".format(run)
    dump = subprocess.Popen([fastq_dump, "--skip-technical", "--stdout", str(archive)], stdout=subprocess.PIPE)
    split = subprocess.Popen(
        [fastp, "--stdin", "--out1", str(staged), "--disable_adapter_trimming", "--disable_quality_filtering",
         "--disable_length_filtering", "--disable_trim_poly_g", "--compression", "4", "--thread", str(threads),
         "--json", str(report), "--html", "/dev/null"],
        stdin=dump.stdout,
    )
    dump.stdout.close()
    split_code, dump_code = split.wait(), dump.wait()
    try:
        if dump_code != 0 or split_code != 0:
            raise RuntimeError("{}: fastq-dump exited {}, fastp exited {}".format(run, dump_code, split_code))
        with staged.open("rb") as handle:
            if handle.read(2) != b"\x1f\x8b":
                raise RuntimeError("{} is not gzip".format(staged))
        total = json.loads(report.read_text())["summary"]["before_filtering"]["total_reads"]
        if total != int(row["read_count"]):
            raise RuntimeError("{}: wrote {} reads, SRA reports {} spots".format(run, total, row["read_count"]))
    except Exception:
        discard([staged, report])
        raise
    report.unlink()
    staged.rename(target)
    print("converted {}".format(run), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locations", type=Path, default=DEFAULT_LOCATIONS)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--sra-dir", type=Path, help="defaults to <analysis-root>/sra")
    parser.add_argument("--verify", action="store_true", help="check the .sra archives against NCBI MD5s")
    parser.add_argument("--convert", action="store_true", help="stream each .sra to gzipped FASTQ")
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
    print("runs: 18 single-end day-41 gill libraries (6 control, 6 high leachate, 6 high tire particles)")
    if args.verify:
        failures = verify_sra(rows, locations, sra_dir, args.workers)
        if failures:
            print("SRA verification failed for: {}".format(sorted(failures)), file=sys.stderr)
            return 2
        print("all {} SRA archives match the NCBI size and MD5".format(len(rows)))
    if args.convert:
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            convert_single_end(row, sra_dir, output_dir, args.threads, args.fastq_dump, args.fastp)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv")
    for test in CONTRASTS:
        write_design_samplesheet([row for row in rows if row["condition"] in ("control", test)],
                                 args.analysis_root / "deseq2_samplesheet_{}.csv".format(test))
    return 0


if __name__ == "__main__":
    sys.exit(main())
