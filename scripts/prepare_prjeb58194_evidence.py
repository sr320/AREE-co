"""Build current-reference processed evidence and mappings for the three PRJEB58194 (PESTO) records."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence
from scripts.prepare_prjna913164_evidence import merge_annotations


ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"
COMMON_FLAGS = ["raw_reanalysis", "pollutant_exposure_response", "pesticide_mixture", "pooled_libraries",
                "stage_blocked_contrast"]
# The two F1 records test against the same six TT libraries.
F1_FLAGS = ["shared_control_group", "dependence_group=PRJEB58194_F1_TT"]
CONTRASTS = {
    "f0_exposed": {
        "study_id": "CGIG_PESTO_F0_RNASEQ_PRJEB58194",
        "sample_comparison": "pesticide_mixture_0_48hpf_vs_control_F0_gastrula_dlarva",
        "quality_flags": [*COMMON_FLAGS, "direct_exposure"],
    },
    "f1_direct": {
        "study_id": "CGIG_PESTO_F1DIRECT_RNASEQ_PRJEB58194",
        "sample_comparison": "pesticide_mixture_TE_vs_TT_F1_gastrula_pediveliger",
        "quality_flags": [*COMMON_FLAGS, *F1_FLAGS, "direct_exposure"],
    },
    "f1_parental": {
        "study_id": "CGIG_PESTO_F1PARENTAL_RNASEQ_PRJEB58194",
        "sample_comparison": "parental_pesticide_mixture_ET_vs_TT_F1_gastrula_pediveliger",
        "quality_flags": [*COMMON_FLAGS, *F1_FLAGS, "intergenerational_exposure"],
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
