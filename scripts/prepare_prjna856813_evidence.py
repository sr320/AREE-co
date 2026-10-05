"""Build current-reference processed evidence and mappings for the two PRJNA856813 tire records."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence
from scripts.prepare_prjna913164_evidence import merge_annotations


ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"
# Both records test against the same six day-41 gill controls.
COMMON_FLAGS = ["raw_reanalysis", "pollutant_exposure_response", "shared_control_group",
                "dependence_group=PRJNA856813_gill_controls"]
CONTRASTS = {
    "leachate_high": {
        "study_id": "CGIG_TIRELEACHATE_RNASEQ_PRJNA856813",
        "sample_comparison": "tire_leachate_high_vs_control_gill_day41",
        "quality_flags": [*COMMON_FLAGS, "tire_leachate"],
    },
    "tire_particles_high": {
        "study_id": "CGIG_TIREPARTICLE_RNASEQ_PRJNA856813",
        "sample_comparison": "tire_microparticles_high_vs_control_gill_day41",
        "quality_flags": [*COMMON_FLAGS, "tire_microparticles"],
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
