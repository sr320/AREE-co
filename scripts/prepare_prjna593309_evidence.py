"""Build current-reference processed evidence and mappings for the two PRJNA593309 records."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"
# The two records use different oysters but the same three 21 C tanks.
COMMON_FLAGS = ["raw_reanalysis", "one_oyster_per_tank", "dependence_group=PRJNA593309_21C_tanks"]
CONTRASTS = {
    "oshv1": {
        "study_id": "CGIG_OSHV1_RNASEQ_PRJNA593309",
        "sample_comparison": "oshv1_cohabitation_48h_vs_0h_21C",
        "quality_flags": ["raw_reanalysis", "pathogen_challenge_response", "unpaired_primary_tank_paired_sensitivity",
                          *COMMON_FLAGS[1:]],
    },
    "heat": {
        "study_id": "CGIG_HEAT_RNASEQ_PRJNA593309",
        "sample_comparison": "29C_vs_21C_12h_oshv1_cohabitation",
        "quality_flags": ["raw_reanalysis", "temperature_response", "during_oshv1_infection", *COMMON_FLAGS[1:]],
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
