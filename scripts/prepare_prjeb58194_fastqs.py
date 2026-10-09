#!/usr/bin/env python3
"""Validate and stage the 54 PRJEB58194 (PESTO) paired-end libraries used in three AREE records.

Sol Dourdin et al. exposed oyster embryos to an environmentally relevant mixture of 18 pesticides
and metals (nominal 2.85 ug/L) from 0 to 48 h post fertilization in two generations. AREE uses:
F0 exposed against control gastrulae and D larvae (blocking on stage); F1 offspring exposed (TE)
against controls (TT); and F1 control offspring of exposed parents (ET) against TT, both blocking
on stage (gastrula, metamorphosis-competent pediveliger). The F0 pediveligers and the F1 EE group
are not used.

All runs come from the NCBI Open Data S3 mirror as ``.sra`` archives listed in
``*_sra_locations.tsv`` (11 of them have no ENA FASTQ): verified against NCBI's published size and
MD5 (``--verify``), then streamed to gzipped FASTQ pairs whose read counts must match the NCBI spot
counts in the manifest (``--convert``), using the steps shared with ``prepare_prjna1329250_fastqs``.

Because it imports the shared SRA steps from ``scripts``, run it from the repository root
as ``python -m scripts.prepare_prjeb58194_fastqs``.
"""

import argparse
import csv
import sys
from pathlib import Path

from aree.raw.fastq import load_run_manifest, sample_name, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjna1329250_fastqs import check_free_space, convert_run, load_locations, verify_sra


DEFAULT_MANIFEST = Path("data/manifests/CGIG_PESTO_RNASEQ_PRJEB58194_runs.tsv")
DEFAULT_LOCATIONS = Path("data/manifests/CGIG_PESTO_RNASEQ_PRJEB58194_sra_locations.tsv")
EXPECTED_CONDITIONS = {"f0_control": 18, "f0_exposed": 18, "f1_TT": 6, "f1_TE": 6, "f1_ET": 6}
# Each record relabels two manifest groups as control and test in its own design sheet.
CONTRASTS = {
    "f0_exposed": ("f0_control", "f0_exposed"),
    "f1_direct": ("f1_TT", "f1_TE"),
    "f1_parental": ("f1_TT", "f1_ET"),
}
# Libraries left out of the DESeq2 design sheets after QC (none yet).
EXCLUDED_SAMPLES = set()


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


def write_contrast_sheet(rows, contrast, path):
    """One record's two groups as control and test; samples keep their manifest names for Salmon."""
    control, test = CONTRASTS[contrast]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "condition", "replicate", "run_accession", "stage", "f0_group"])
        writer.writeheader()
        for row in rows:
            if row["condition"] in (control, test) and sample_name(row) not in EXCLUDED_SAMPLES:
                writer.writerow({"sample": sample_name(row), "condition": "control" if row["condition"] == control else contrast,
                                 "replicate": row["replicate"], "run_accession": row["run_accession"],
                                 "stage": row["stage"], "f0_group": row["f0_group"]})
    print("wrote DESeq2 design sheet {}".format(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locations", type=Path, default=DEFAULT_LOCATIONS)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--sra-dir", type=Path, help="defaults to <analysis-root>/sra")
    parser.add_argument("--verify", action="store_true", help="check the .sra archives against NCBI MD5s")
    parser.add_argument("--convert", action="store_true", help="stream each .sra to a gzipped FASTQ pair")
    parser.add_argument("--keep-sra", action="store_true", help="retain .sra archives after conversion")
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
    print("runs: 54 paired-end embryo and larva libraries (F0 control/exposed x 2 stages x 9; F1 TT/TE/ET x 2 stages x 3)")
    if args.verify:
        failures = verify_sra(rows, locations, sra_dir, args.workers)
        if failures:
            print("SRA verification failed for: {}".format(sorted(failures)), file=sys.stderr)
            return 2
        print("all {} SRA archives match the NCBI size and MD5".format(len(rows)))
    if args.convert:
        shortfall = check_free_space(rows, sra_dir, output_dir, args.keep_sra)
        if shortfall:
            print(shortfall, file=sys.stderr)
            return 2
        for row in sorted(rows, key=lambda item: item["run_accession"]):
            convert_run(row, sra_dir, output_dir, args.threads, args.fastq_dump, args.fastp, args.keep_sra)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("stage", "generation", "f0_group"))
    for test in CONTRASTS:
        write_contrast_sheet(rows, test, args.analysis_root / "deseq2_samplesheet_{}.csv".format(test))
    return 0


if __name__ == "__main__":
    sys.exit(main())
