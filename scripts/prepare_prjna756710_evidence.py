"""Build current-reference processed evidence and mappings for PRJNA756710."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


SAMPLE_COMPARISON = "low_salinity_10ppt_vs_30ppt_during_vibrio_infection"
QUALITY_FLAGS = ["raw_reanalysis", "low_salinity_response", "infection_background_all_arms", "timepoint_blocked_contrast"]
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
