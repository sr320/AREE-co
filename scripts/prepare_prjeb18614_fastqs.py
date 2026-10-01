#!/usr/bin/env python3
"""Validate and stage the 26 PRJEB18614 gill FASTQ pairs (one lane per library).

Reads come from two mirrors. Runs listed in ``*_sra_locations.tsv`` are taken from the NCBI
Open Data S3 mirror as ``.sra`` archives: verified against NCBI's published size and MD5
(``--verify``), then streamed to gzipped FASTQ pairs whose read counts must match the
manifest (``--convert``), using the steps shared with ``prepare_prjna1329250_fastqs``.
ENA slowed this transfer to about 0.5 MB/s and refused requests with HTTP 403 partway
through, so the runs it had not delivered intact were moved to NCBI.

Every other run comes from the ENA FASTQ mirror and must match the manifest size and MD5.
``--download`` fetches those with a single resumable connection. For a faster transfer,
fetch the manifest URLs with a multi-connection client into ``<analysis-root>/fastq`` as
``<run>_<mate>.fastq.gz.part`` and then run ``--download``: a complete ``.part`` file is
promoted only after it matches the manifest size and MD5. Use one connection per file:
with four, about 60% of these files arrived corrupted.

Because it imports the shared SRA steps from ``scripts``, run it from the repository root
as ``python -m scripts.prepare_prjeb18614_fastqs``.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import download_runs, load_run_manifest, preflight, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_HAB_RNASEQ_PRJEB18614_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_HAB_RNASEQ_PRJEB18614_sra_locations.tsv")
EXPECTED_CONDITIONS = {"control": 13, "alexandrium": 13}


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def split_by_source(rows, locations):
    """Split manifest rows into (ENA-sourced, NCBI-sourced) by the SRA locations table."""
    unknown = sorted(set(locations) - {row["run_accession"] for row in rows})
    if unknown:
        raise ValueError("SRA locations name runs that are not in the manifest: {}".format(unknown))
    ena = [row for row in rows if row["run_accession"] not in locations]
    ncbi = [row for row in rows if row["run_accession"] in locations]
    return ena, ncbi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locations", type=Path, default=DEFAULT_LOCATIONS)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--download", action="store_true", help="fetch or promote the ENA-sourced runs")
    parser.add_argument("--sra-dir", type=Path, help="defaults to <analysis-root>/sra")
    parser.add_argument("--verify", action="store_true", help="check the NCBI .sra archives against NCBI MD5s")
    parser.add_argument("--convert", action="store_true", help="stream each NCBI .sra to a gzipped FASTQ pair")
    parser.add_argument("--keep-sra", action="store_true", help="retain .sra archives after conversion")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--retries", type=int, default=12)
    parser.add_argument("--threads", type=int, default=8, help="fastp compression threads")
    parser.add_argument("--fastq-dump", default="fastq-dump")
    parser.add_argument("--fastp", default="fastp")
    parser.add_argument("--scratch-multiplier", type=float, default=1.25)
    args = parser.parse_args()
    rows = load_manifest(args.manifest)
    locations = load_locations(args.locations)
    ena, ncbi = split_by_source(rows, locations)
    output_dir = args.analysis_root / "fastq"
    output_dir.mkdir(parents=True, exist_ok=True)
    sra_dir = args.sra_dir or args.analysis_root / "sra"
    print("runs: 26 (13 H. triquetra controls, 13 A. minutum; one pool per diet at each of 13 times)")
    print("sources: {} from ENA FASTQs, {} from NCBI .sra archives".format(len(ena), len(ncbi)))
    if not preflight(ena, args.analysis_root, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    if args.download:
        download_runs(ena, output_dir, workers=args.workers, retries=args.retries, timeout=90)
    if args.verify:
        failures = verify_sra(ncbi, locations, sra_dir, args.workers)
        if failures:
            print("SRA verification failed for: {}".format(sorted(failures)), file=sys.stderr)
            return 2
        print("all {} SRA archives match the NCBI size and MD5".format(len(ncbi)))
    if args.convert:
        shortfall = check_free_space(ncbi, sra_dir, output_dir, args.keep_sra)
        if shortfall:
            print(shortfall, file=sys.stderr)
            return 2
        for row in sorted(ncbi, key=lambda item: item["run_accession"]):
            convert_run(row, sra_dir, output_dir, args.threads, args.fastq_dump, args.fastp, args.keep_sra)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    # timepoint pairs the diets; sequencing_pool feeds the multiplex-pool sensitivity model.
    write_design_samplesheet(
        rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("timepoint", "sequencing_pool"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
