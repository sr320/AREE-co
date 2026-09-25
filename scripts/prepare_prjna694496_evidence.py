"""Build processed evidence, exact-ID mappings, and RefSeq annotations for PRJNA694496."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path
from urllib.parse import unquote


PROCESSED_FIELDS = [
    "sample_comparison", "feature_id_original", "feature_type", "molecular_direction",
    "effect_size", "effect_size_type", "standard_error", "p_value", "adjusted_p_value",
    "analysis_method", "quality_flags",
]
MAPPING_FIELDS = [
    "feature_id_original", "feature_id_standardized", "ortholog_reference",
    "mapping_confidence", "mapping_release", "mapping_evidence",
]
ANNOTATION_FIELDS = [
    "feature_id_standardized", "gene_symbol", "description", "gene_biotype",
    "annotation_release",
]


def parse_attributes(text):
    attributes = {}
    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attributes[key] = unquote(value)
    return attributes


def load_gff_annotations(path):
    annotations = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "pseudogene"}:
                continue
            attributes = parse_attributes(fields[8])
            match = re.search(r"(?:^|,)GeneID:(\d+)(?:,|$)", attributes.get("Dbxref", ""))
            if not match:
                continue
            standardized = "NCBI:GeneID:" + match.group(1)
            annotations[standardized] = {
                "feature_id_standardized": standardized,
                "gene_symbol": attributes.get("gene", attributes.get("Name", "")),
                "description": attributes.get("description", ""),
                "gene_biotype": attributes.get("gene_biotype", ""),
                "annotation_release": "GCF_963853765.1-RS_2024_06",
            }
    return annotations


def prepare(results_path, gff_path, processed_path, mapping_path, annotation_path):
    results = list(csv.DictReader(results_path.open(newline=""), delimiter="\t"))
    required = {"feature_id_standardized", "log2FoldChange", "lfcSE", "pvalue", "padj"}
    if not results or not required.issubset(results[0]):
        raise ValueError("Gene results are empty or missing required DESeq2 columns")
    if len({row["feature_id_standardized"] for row in results}) != len(results):
        raise ValueError("Gene results contain duplicate standardized identifiers")
    annotations = load_gff_annotations(gff_path)

    processed_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    annotation_path.parent.mkdir(parents=True, exist_ok=True)
    with processed_path.open("w", newline="") as processed_handle, \
            mapping_path.open("w", newline="") as mapping_handle, \
            annotation_path.open("w", newline="") as annotation_handle:
        processed_writer = csv.DictWriter(processed_handle, fieldnames=PROCESSED_FIELDS, delimiter="\t")
        mapping_writer = csv.DictWriter(mapping_handle, fieldnames=MAPPING_FIELDS, delimiter="\t")
        annotation_writer = csv.DictWriter(annotation_handle, fieldnames=ANNOTATION_FIELDS, delimiter="\t")
        processed_writer.writeheader()
        mapping_writer.writeheader()
        annotation_writer.writeheader()
        missing_annotations = []
        for row in results:
            identifier = row["feature_id_standardized"]
            if not re.fullmatch(r"NCBI:GeneID:\d+", identifier):
                raise ValueError(f"Unexpected current-reference identifier: {identifier}")
            effect = float(row["log2FoldChange"])
            flags = [
                "raw_reanalysis", "qc_reviewed", "population_selection_association",
                "constitutive_expression", "selected_A_pc2_displacement",
                "sensitivity_direction_robust",
            ]
            if row["padj"] and float(row["padj"]) < 0.05:
                flags.append("fdr_lt_0.05")
            if not row["pvalue"] or not row["padj"]:
                flags.append("deseq2_significance_unavailable")
            processed_writer.writerow({
                "sample_comparison": "thermotolerance_selected_vs_nonselected_control",
                "feature_id_original": identifier,
                "feature_type": "gene",
                "molecular_direction": "up" if effect > 0 else "down" if effect < 0 else "unknown",
                "effect_size": row["log2FoldChange"],
                "effect_size_type": "log2_fold_change",
                "standard_error": row["lfcSE"],
                "p_value": row["pvalue"],
                "adjusted_p_value": row["padj"],
                "analysis_method": "Salmon_1.10.3_tximport_1.30.0_DESeq2_1.42.0_unshrunk_effect",
                "quality_flags": ";".join(flags),
            })
            mapping_writer.writerow({
                "feature_id_original": identifier,
                "feature_id_standardized": identifier,
                "ortholog_reference": "",
                "mapping_confidence": "exact",
                "mapping_release": "GCF_963853765.1-RS_2024_06",
                "mapping_evidence": "direct_gene_id_from_versioned_refseq_tx2gene",
            })
            annotation = annotations.get(identifier)
            if annotation is None:
                missing_annotations.append(identifier)
                annotation = {
                    "feature_id_standardized": identifier, "gene_symbol": "",
                    "description": "", "gene_biotype": "",
                    "annotation_release": "GCF_963853765.1-RS_2024_06",
                }
            annotation_writer.writerow(annotation)
    if missing_annotations:
        raise ValueError(
            f"{len(missing_annotations)} tested GeneIDs are absent from the versioned GFF; "
            f"first: {missing_annotations[0]}"
        )
    print(f"wrote {len(results)} processed effects, exact mappings, and RefSeq annotations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--processed", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.results, args.gff, args.processed, args.mapping, args.annotations)
