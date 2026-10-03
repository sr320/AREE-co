"""Build current-reference processed evidence and mappings for PRJNA735889 (Lutier 2022 juvenile pH gradient)."""

import argparse
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


SAMPLE_COMPARISON = "ph_6.5_6.8_vs_ph_7.4_7.8_23d_juveniles"
QUALITY_FLAGS = ["raw_reanalysis", "ocean_acidification_response", "ph_band_contrast_tipping_window_excluded", "no_ambient_ph_reference", "one_tank_per_ph_level", "gametogenesis_blocked_contrast"]
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
