"""Export DESeq2 gene results as AREE processed evidence, exact-ID mappings and RefSeq annotations."""

import csv
import re
from pathlib import Path

from aree.raw.gff import gene_annotations


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
REQUIRED_RESULT_COLUMNS = {"feature_id_standardized", "log2FoldChange", "lfcSE", "pvalue", "padj"}
REFSEQ_RELEASE = "GCF_963853765.1-RS_2024_06"
ANALYSIS_METHOD = "Salmon_1.10.3_tximport_1.30.0_DESeq2_1.42.0_unshrunk_effect"


def _clean_numeric(value):
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text == "" or text.lower() in {"na", "nan"} else text


def export_deseq2_evidence(
    results_path,
    gff_path,
    processed_path,
    mapping_path,
    annotation_path,
    sample_comparison,
    quality_flags,
    reference_release=REFSEQ_RELEASE,
    analysis_method=ANALYSIS_METHOD,
):
    """Write the three AREE input tables for one DESeq2 contrast on a versioned RefSeq reference.

    ``quality_flags`` are the study-level flags; ``fdr_lt_0.05`` and
    ``deseq2_significance_unavailable`` are added per gene. Every tested GeneID must be present
    in the GFF; otherwise a ValueError names the first missing one.
    """
    with results_path.open(newline="") as handle:
        results = list(csv.DictReader(handle, delimiter="\t"))
    if not results or not REQUIRED_RESULT_COLUMNS.issubset(results[0]):
        raise ValueError("Gene results are empty or missing required DESeq2 columns")
    if len({row["feature_id_standardized"] for row in results}) != len(results):
        raise ValueError("Gene results contain duplicate standardized identifiers")
    annotations = gene_annotations(gff_path, reference_release)
    for path in (processed_path, mapping_path, annotation_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    missing = []
    with processed_path.open("w", newline="") as processed_handle, \
            mapping_path.open("w", newline="") as mapping_handle, \
            annotation_path.open("w", newline="") as annotation_handle:
        processed_writer = csv.DictWriter(processed_handle, fieldnames=PROCESSED_FIELDS, delimiter="\t")
        mapping_writer = csv.DictWriter(mapping_handle, fieldnames=MAPPING_FIELDS, delimiter="\t")
        annotation_writer = csv.DictWriter(annotation_handle, fieldnames=ANNOTATION_FIELDS, delimiter="\t")
        processed_writer.writeheader()
        mapping_writer.writeheader()
        annotation_writer.writeheader()
        for row in results:
            identifier = row["feature_id_standardized"]
            if not re.fullmatch(r"NCBI:GeneID:\d+", identifier):
                raise ValueError("Unexpected current-reference identifier: {}".format(identifier))
            effect = float(row["log2FoldChange"])
            se = _clean_numeric(row["lfcSE"])
            pvalue = _clean_numeric(row["pvalue"])
            padj = _clean_numeric(row["padj"])
            flags = list(quality_flags)
            if padj and float(padj) < 0.05:
                flags.append("fdr_lt_0.05")
            if not pvalue or not padj:
                flags.append("deseq2_significance_unavailable")
            processed_writer.writerow({
                "sample_comparison": sample_comparison,
                "feature_id_original": identifier,
                "feature_type": "gene",
                "molecular_direction": "up" if effect > 0 else "down" if effect < 0 else "unknown",
                "effect_size": row["log2FoldChange"],
                "effect_size_type": "log2_fold_change",
                "standard_error": se,
                "p_value": pvalue,
                "adjusted_p_value": padj,
                "analysis_method": analysis_method,
                "quality_flags": ";".join(flags),
            })
            mapping_writer.writerow({
                "feature_id_original": identifier,
                "feature_id_standardized": identifier,
                "ortholog_reference": "",
                "mapping_confidence": "exact",
                "mapping_release": reference_release,
                "mapping_evidence": "direct_gene_id_from_versioned_refseq_tx2gene",
            })
            annotation = annotations.get(identifier)
            if annotation is None:
                missing.append(identifier)
            else:
                annotation_writer.writerow(annotation)
    if missing:
        raise ValueError(
            "{} tested GeneIDs are absent from the versioned GFF; first: {}".format(len(missing), missing[0])
        )
    print("wrote {} processed effects, exact mappings, and RefSeq annotations".format(len(results)))
    return len(results)


def merge_annotations(paths, output):
    """Union the per-contrast RefSeq annotation tables into the one table keyed by bioproject."""
    header, rows = None, {}
    for path in paths:
        with Path(path).open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            current = next(reader)
            if header is None:
                header = current
            elif current != header:
                raise ValueError("annotation tables have different columns: {}".format(path))
            key = header.index("feature_id_standardized")
            for row in reader:
                rows.setdefault(row[key], row)
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows[feature] for feature in sorted(rows))
