"""Build current-reference processed evidence and mappings for PRJNA516762."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


SAMPLE_COMPARISON = "heat_35C_2h_vs_control_12C"
QUALITY_FLAGS = ["raw_reanalysis", "acute_heat_response", "pooled_libraries"]


def prepare(results_path, gff_path, processed_path, mapping_path, annotation_path):
    export_deseq2_evidence(
        results_path, gff_path, processed_path, mapping_path, annotation_path,
        sample_comparison=SAMPLE_COMPARISON, quality_flags=QUALITY_FLAGS,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.results, args.gff, args.processed, args.mapping, args.annotations)
