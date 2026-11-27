"""Build current-reference processed evidence and mappings for PRJNA735889 (Lutier 2022 juvenile pH gradient)."""

import argparse
import csv
import gzip
import tempfile
from pathlib import Path

from aree.raw.deseq2_evidence import export_deseq2_evidence


SAMPLE_COMPARISON = "ph_6.5_6.8_vs_ph_7.4_7.8_23d_juveniles"
QUALITY_FLAGS = [
    "raw_reanalysis",
    "ocean_acidification_response",
    "ph_band_contrast_tipping_window_excluded",
    "no_ambient_ph_reference",
    "one_tank_per_ph_level",
    "gametogenesis_blocked_contrast",
    "tank_random_effect_mixed_model",
]
ANALYSIS_METHOD = "Salmon_2.3.4_tximport_1.30.0_limma_3.58.1_dream_1.32.5_tank_random_effect"


def prepare(results_path, gff_path, processed_path, mapping_path, annotation_path, existing_annotations=None):
    """Export PRJNA735889 results to processed evidence, exact mapping, and RefSeq annotations."""
    if gff_path is not None:
        export_deseq2_evidence(
            results_path, gff_path, processed_path, mapping_path, annotation_path,
            sample_comparison=SAMPLE_COMPARISON, quality_flags=QUALITY_FLAGS,
            analysis_method=ANALYSIS_METHOD,
        )
        return
    if existing_annotations is None:
        raise ValueError("Provide either --gff or --existing-annotations")
    # Build a minimal GFF from the versioned RefSeq annotation table so export_deseq2_evidence
    # performs identical GeneID validation and writes all three outputs consistently.
    with existing_annotations.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with tempfile.NamedTemporaryFile(suffix=".gff.gz", delete=False) as tmp:
        tmp_gff = Path(tmp.name)
    try:
        with gzip.open(tmp_gff, "wt") as gz:
            for row in rows:
                gid = row["feature_id_standardized"].split(":")[-1]
                desc = row["description"].replace("%", "%25").replace(";", "%3B").replace("=", "%3D")
                sym = row["gene_symbol"].replace("%", "%25").replace(";", "%3B").replace("=", "%3D")
                bio = row["gene_biotype"]
                gz.write(
                    "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
                    "ID=gene-{0};Dbxref=GeneID:{0};Name={1};gene={1};description={2};gene_biotype={3}\n".format(
                        gid, sym, desc, bio
                    )
                )
        export_deseq2_evidence(
            results_path, tmp_gff, processed_path, mapping_path, annotation_path,
            sample_comparison=SAMPLE_COMPARISON, quality_flags=QUALITY_FLAGS,
            analysis_method=ANALYSIS_METHOD,
        )
    finally:
        if tmp_gff.exists():
            tmp_gff.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--gff", type=Path, default=None)
    parser.add_argument("--existing-annotations", type=Path, default=None)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    args = parser.parse_args()
    prepare(
        args.results, args.gff, args.processed, args.mapping, args.annotations,
        existing_annotations=args.existing_annotations,
    )
