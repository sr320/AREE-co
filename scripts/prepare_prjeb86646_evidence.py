"""Build current-reference processed evidence and mappings for PRJEB86646 (DECICOMP 30 C)."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


SAMPLE_COMPARISON = "30C_vs_23C_fed_before_oshv1"
QUALITY_FLAGS = ["raw_reanalysis", "temperature_response", "pre_infection_baseline", "family_blocked_contrast", "dependence_group=DECICOMP_PRJEB86525_0h_controls", "shared_control_group_with_PRJEB86618"]
ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"


def prepare(results_path, gff_path, processed_path, mapping_path, annotation_path):
    export_deseq2_evidence(
        results_path, gff_path, processed_path, mapping_path, annotation_path,
        sample_comparison=SAMPLE_COMPARISON, quality_flags=QUALITY_FLAGS,
        analysis_method=ANALYSIS_METHOD,
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
