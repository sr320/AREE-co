"""Build current-reference processed evidence and mappings for the three PRJNA298285 larval records."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence
from scripts.prepare_prjna913164_evidence import merge_annotations


ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"
# All three records test against the same six ambient larval libraries.
COMMON_FLAGS = ["raw_reanalysis", "larval_development", "pooled_libraries", "stage_blocked_contrast",
                "shared_control_group", "dependence_group=PRJNA298285_ambient_larvae", "no_linked_publication"]
CONTRASTS = {
    "acidified": {
        "study_id": "CGIG_LARVAL_OA_RNASEQ_PRJNA298285",
        "sample_comparison": "pH7.9_vs_pH8.1_20C_larvae",
        "quality_flags": [*COMMON_FLAGS, "ocean_acidification_response"],
    },
    "warmed": {
        "study_id": "CGIG_LARVAL_HEAT_RNASEQ_PRJNA298285",
        "sample_comparison": "22C_vs_20C_pH8.1_larvae",
        "quality_flags": [*COMMON_FLAGS, "temperature_response"],
    },
    "acidified_warmed": {
        "study_id": "CGIG_LARVAL_MULTISTRESS_RNASEQ_PRJNA298285",
        "sample_comparison": "pH7.9_22C_vs_pH8.1_20C_larvae",
        "quality_flags": [*COMMON_FLAGS, "multiple_stressor_response"],
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
