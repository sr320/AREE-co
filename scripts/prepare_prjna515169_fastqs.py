#!/usr/bin/env python3
"""Validate and stage the 9 PRJNA515169 paired-end libraries used in two AREE Vibrio records.

Rubio et al. 2019 ran three independent infection experiments. In each, 30 spat per condition were
injected with a Vibrio strain or with sterile seawater (mock), and one pool of 30 whole oysters per
condition was sequenced 8 h after injection. AREE tests each virulent strain, V. crassostreae J2-9
and V. tasmaniensis LGP32, against the shared mock pools, blocking on experiment. The nonvirulent
strains and the untreated oysters are not used.

``--download`` fetches each file with a single resumable connection. For a faster transfer,
fetch the manifest URLs with a multi-connection client into ``<analysis-root>/fastq`` as
``<run>_<mate>.fastq.gz.part`` and then run ``--download``: a complete ``.part`` file is
promoted only after it matches the manifest size and MD5. Use one connection per file.
"""

import argparse
import sys
from pathlib import Path

from aree.raw.fastq import download_runs, load_run_manifest, preflight, write_design_samplesheet, write_nfcore_samplesheet


DEFAULT_MANIFEST = Path("data/manifests/CGIG_VIBRIO_RNASEQ_PRJNA515169_runs.tsv")
EXPECTED_CONDITIONS = {"control": 3, "vcrassostreae_J2_9": 3, "vtasmaniensis_LGP32": 3}
CONTRASTS = ("vcrassostreae_J2_9", "vtasmaniensis_LGP32")


def load_manifest(path):
    return load_run_manifest(path, EXPECTED_CONDITIONS)


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
    print("runs: 9 paired-end pools of 30 spat (mock, V. crassostreae J2-9, V. tasmaniensis LGP32; 3 experiments)")
    if not preflight(rows, args.analysis_root, args.scratch_multiplier):
        print("preflight failed: choose a scratch location with more free space", file=sys.stderr)
        return 2
    if args.download:
        download_runs(rows, output_dir, workers=args.workers, retries=args.retries, timeout=90)
    write_nfcore_samplesheet(rows, output_dir, args.analysis_root / "nfcore_samplesheet.csv")
    write_design_samplesheet(rows, args.analysis_root / "deseq2_samplesheet.csv", extra_columns=("experiment",))
    for test in CONTRASTS:
        write_design_samplesheet([row for row in rows if row["condition"] in ("control", test)],
                                 args.analysis_root / "deseq2_samplesheet_{}.csv".format(test), extra_columns=("experiment",))
    return 0


if __name__ == "__main__":
    sys.exit(main())
