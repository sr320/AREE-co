"""Build current-reference processed evidence and mappings for the two PRJNA407831 heat records."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence
from scripts.prepare_prjna913164_evidence import merge_annotations


ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"
# Both records test against the same twelve 0 h pools.
COMMON_FLAGS = ["raw_reanalysis", "acute_heat_response", "pooled_libraries", "population_blocked_contrast",
                "shared_control_group", "dependence_group=PRJNA407831_0h_pools"]
CONTRASTS = {
    "h6": {
        "study_id": "CGIG_HEAT6H_RNASEQ_PRJNA407831",
        "sample_comparison": "35C_6h_vs_0h_gill",
        "quality_flags": [*COMMON_FLAGS, "heat_6h"],
    },
    "h24": {
        "study_id": "CGIG_HEAT24H_RNASEQ_PRJNA407831",
        "sample_comparison": "35C_24h_vs_0h_gill",
        "quality_flags": [*COMMON_FLAGS, "heat_24h"],
    },
}


def prepare(contrast, results_path, gff_path, processed_path, mapping_path, annotation_path):
    spec = CONTRASTS[contrast]
    export_deseq2_evidence(
        results_path, gff_path, processed_path, mapping_path, annotation_path,
        sample_comparison=spec["sample_comparison"], quality_flags=spec["quality_flags"],
        analysis_method=ANALYSIS_METHOD,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contrast", choices=sorted(CONTRASTS))
    parser.add_argument("--results", type=Path)
    parser.add_argument("--gff", type=Path)
    parser.add_argument("--processed", type=Path)
    parser.add_argument("--mapping", type=Path)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--merge", type=Path, nargs="+",
                        help="instead of exporting, union these annotation tables into --annotations")
    args = parser.parse_args()
    if args.merge:
        merge_annotations(args.merge, args.annotations)
    else:
        prepare(args.contrast, args.results, args.gff, args.processed, args.mapping, args.annotations)
