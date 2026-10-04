"""Build current-reference processed evidence and mappings for the two PRJNA913164 heatwave records."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


ANALYSIS_METHOD = "Salmon_2.3.4_noLengthCorrection_tximport_1.30.0_DESeq2_1.42.1_tagseq_counts_unshrunk_effect"
COMMON_FLAGS = ["raw_reanalysis", "quantseq_3prime_tag_seq", "ploidy_blocked_contrast", "shared_control_group"]
# Each stressor is a separate study record tested against the shared 20 C controls.
CONTRASTS = {
    "heat": {
        "study_id": "CGIG_HEAT_RNASEQ_PRJNA913164",
        "sample_comparison": "30C_seawater_vs_20C_adult_ctenidium",
        "quality_flags": [*COMMON_FLAGS, "temperature_response"],
    },
    "heat_emersion": {
        "study_id": "CGIG_MULTISTRESS_RNASEQ_PRJNA913164",
        "sample_comparison": "30C_seawater_then_44C_emersion_6h_vs_20C_adult_ctenidium",
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
    parser.add_argument("--contrast", choices=sorted(CONTRASTS), required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.contrast, args.results, args.gff, args.processed, args.mapping, args.annotations)
