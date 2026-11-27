"""Validate and export multi-assay differential/quantitative tables to AREE processed evidence."""

import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from aree.harmonize.identifiers import CONFIDENCE_ORDER
from aree.harmonize.processed import REQUIRED_PROCESSED_COLUMNS
from aree.raw.deseq2_evidence import (
    ANNOTATION_FIELDS,
    MAPPING_FIELDS,
    PROCESSED_FIELDS,
    _clean_numeric,
)

ALLOWED_DIRECTIONS = {
    "up",
    "down",
    "higher",
    "lower",
    "hypermethylated",
    "hypomethylated",
    "mixed",
    "unknown",
}

ALLOWED_CONFIDENCES = set(CONFIDENCE_ORDER)

DIRECTION_LABELS = {
    "expression": ("up", "down"),
    "methylation": ("hypermethylated", "hypomethylated"),
    "abundance": ("higher", "lower"),
}


def _parse_optional_float(value: str, column: str, row_idx: int) -> Optional[float]:
    text = "" if value is None else str(value).strip()
    if text == "" or text.lower() in {"na", "nan"}:
        return None
    try:
        parsed = float(text)
    except ValueError as exc:
        raise ValueError("Row {}: invalid float {!r} in column {}".format(row_idx, value, column)) from exc
    if not math.isfinite(parsed):
        raise ValueError("Row {}: non-finite value {!r} in column {}".format(row_idx, value, column))
    return parsed


def validate_processed_table(
    processed_path: Path,
    mapping_path: Optional[Path] = None,
    annotation_path: Optional[Path] = None,
    out_path: Optional[Path] = None,
) -> dict:
    """Validate exported processed evidence, identifier mappings, and annotations against AREE constraints."""
    processed_path = Path(processed_path)
    with processed_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = list(reader.fieldnames or [])
        rows = list(reader)

    missing_cols = [col for col in REQUIRED_PROCESSED_COLUMNS if col not in columns]
    if missing_cols:
        raise ValueError("Processed evidence missing required columns: {}".format(", ".join(missing_cols)))
    if not rows:
        raise ValueError("Processed evidence table is empty: {}".format(processed_path))

    fdr_significant = 0
    feature_types = set()
    original_ids = set()
    for idx, row in enumerate(rows, start=1):
        for col in REQUIRED_PROCESSED_COLUMNS:
            val = row.get(col)
            if val is None or str(val).strip() == "":
                if col not in {"standard_error", "p_value", "adjusted_p_value"}:
                    raise ValueError("Row {}: missing required value in {}".format(idx, col))
        direction = row["molecular_direction"].strip()
        if direction not in ALLOWED_DIRECTIONS:
            raise ValueError("Row {}: invalid molecular_direction {!r}".format(idx, direction))
        _parse_optional_float(row["effect_size"], "effect_size", idx)
        se = _parse_optional_float(row.get("standard_error", ""), "standard_error", idx)
        if se is not None and se < 0:
            raise ValueError("Row {}: negative standard_error {}".format(idx, se))
        pval = _parse_optional_float(row.get("p_value", ""), "p_value", idx)
        if pval is not None and not (0.0 <= pval <= 1.0):
            raise ValueError("Row {}: p_value {} outside [0, 1]".format(idx, pval))
        padj = _parse_optional_float(row.get("adjusted_p_value", ""), "adjusted_p_value", idx)
        if padj is not None:
            if not (0.0 <= padj <= 1.0):
                raise ValueError("Row {}: adjusted_p_value {} outside [0, 1]".format(idx, padj))
            if padj < 0.05:
                fdr_significant += 1
        feature_types.add(row["feature_type"].strip())
        original_ids.add(row["feature_id_original"].strip())

    mapping_rows = 0
    if mapping_path is not None:
        mapping_path = Path(mapping_path)
        with mapping_path.open(newline="", encoding="utf-8") as handle:
            m_reader = csv.DictReader(handle, delimiter="\t")
            m_cols = set(m_reader.fieldnames or [])
            m_list = list(m_reader)
        required_m = {"feature_id_original", "feature_id_standardized", "mapping_confidence"}
        if not required_m.issubset(m_cols):
            raise ValueError("Mapping table missing columns: {}".format(", ".join(sorted(required_m - m_cols))))
        mapped_orig = set()
        for idx, row in enumerate(m_list, start=1):
            conf = (row.get("mapping_confidence") or "").strip()
            if conf not in ALLOWED_CONFIDENCES:
                raise ValueError("Mapping row {}: invalid mapping_confidence {!r}".format(idx, conf))
            mapped_orig.add((row.get("feature_id_original") or "").strip())
        missing_map = original_ids - mapped_orig
        if missing_map:
            sample_missing = ", ".join(sorted(missing_map)[:5])
            raise ValueError("Processed features missing from mapping table: {}".format(sample_missing))
        mapping_rows = len(m_list)

    annotation_rows = 0
    if annotation_path is not None:
        annotation_path = Path(annotation_path)
        with annotation_path.open(newline="", encoding="utf-8") as handle:
            a_reader = csv.DictReader(handle, delimiter="\t")
            a_cols = set(a_reader.fieldnames or [])
            a_list = list(a_reader)
        required_a = {"feature_id_standardized", "gene_symbol", "description", "gene_biotype", "annotation_release"}
        if not required_a.issubset(a_cols):
            raise ValueError("Annotation table missing columns: {}".format(", ".join(sorted(required_a - a_cols))))
        annotation_rows = len(a_list)

    report = {
        "status": "valid",
        "processed_path": str(processed_path),
        "n_rows": len(rows),
        "n_unique_features": len(original_ids),
        "n_fdr_lt_0_05": fdr_significant,
        "feature_types": sorted(feature_types),
        "mapping_rows": mapping_rows,
        "annotation_rows": annotation_rows,
    }
    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def _first_present(row: Dict[str, str], candidates: List[str], default: str = "") -> str:
    for key in candidates:
        if key in row and str(row[key]).strip() != "":
            return str(row[key]).strip()
    return default


def export_assay_evidence(
    results_path: Path,
    processed_path: Path,
    mapping_path: Path,
    annotation_path: Path,
    sample_comparison: str,
    feature_type: str,
    effect_size_type: str,
    analysis_method: str,
    quality_flags: Optional[Iterable[str]] = None,
    reference_release: str = "assay_feature_catalog_v1",
    direction_mode: str = "abundance",
    ortholog_reference: str = "none",
    mapping_release: str = "direct_assay_id_v1",
    mapping_evidence: str = "Assay-native standardized feature identifier",
) -> dict:
    """Export generic differential/quantitative assay results (WGBS, proteomics, metabolomics) to AREE tables."""
    if direction_mode not in DIRECTION_LABELS:
        raise ValueError("Unsupported direction_mode {!r}".format(direction_mode))
    pos_label, neg_label = DIRECTION_LABELS[direction_mode]

    results_path = Path(results_path)
    with results_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Assay results table is empty: {}".format(results_path))

    study_flags = [f.strip() for f in (quality_flags or []) if f and f.strip()]
    for path in (processed_path, mapping_path, annotation_path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    seen_original = set()
    seen_annotations = set()

    with Path(processed_path).open("w", newline="", encoding="utf-8") as p_handle, \
            Path(mapping_path).open("w", newline="", encoding="utf-8") as m_handle, \
            Path(annotation_path).open("w", newline="", encoding="utf-8") as a_handle:
        p_writer = csv.DictWriter(p_handle, fieldnames=PROCESSED_FIELDS, delimiter="\t")
        m_writer = csv.DictWriter(m_handle, fieldnames=MAPPING_FIELDS, delimiter="\t")
        a_writer = csv.DictWriter(a_handle, fieldnames=ANNOTATION_FIELDS, delimiter="\t")
        p_writer.writeheader()
        m_writer.writeheader()
        a_writer.writeheader()

        for idx, row in enumerate(rows, start=1):
            orig_id = _first_present(row, ["feature_id_original", "feature_id", "id", "site_id", "protein_id", "metabolite_id"])
            if not orig_id:
                raise ValueError("Row {}: missing feature identifier column".format(idx))
            if orig_id in seen_original:
                raise ValueError("Duplicate feature identifier in assay results: {}".format(orig_id))
            seen_original.add(orig_id)

            std_id = _first_present(row, ["feature_id_standardized"], default=orig_id)
            effect_raw = _first_present(row, ["effect_size", "log2FoldChange", "logFC", "delta_beta", "mean_diff"])
            if not effect_raw:
                raise ValueError("Row {}: missing effect size column".format(idx))
            effect = float(effect_raw)
            se_raw = _clean_numeric(_first_present(row, ["standard_error", "lfcSE", "se"]))
            if not se_raw:
                t_val = _clean_numeric(_first_present(row, ["t_stat", "t"]))
                if t_val and float(t_val) != 0.0 and effect != 0.0:
                    se_raw = "{:.6g}".format(abs(effect / float(t_val)))
            pval = _clean_numeric(_first_present(row, ["p_value", "P.Value", "pvalue"]))
            padj = _clean_numeric(_first_present(row, ["adjusted_p_value", "adj.P.Val", "padj", "fdr"]))

            row_flags = list(study_flags)
            if padj and float(padj) < 0.05:
                row_flags.append("fdr_lt_0.05")
            if not pval or not padj:
                row_flags.append("significance_unavailable")

            direction = pos_label if effect > 0 else neg_label if effect < 0 else "unknown"
            p_writer.writerow({
                "sample_comparison": sample_comparison,
                "feature_id_original": orig_id,
                "feature_type": feature_type,
                "molecular_direction": direction,
                "effect_size": effect_raw,
                "effect_size_type": effect_size_type,
                "standard_error": se_raw,
                "p_value": pval,
                "adjusted_p_value": padj,
                "analysis_method": analysis_method,
                "quality_flags": ";".join(row_flags),
            })

            conf = _first_present(row, ["mapping_confidence"], default="exact" if std_id == orig_id else "inferred")
            m_writer.writerow({
                "feature_id_original": orig_id,
                "feature_id_standardized": std_id,
                "ortholog_reference": _first_present(row, ["ortholog_reference"], default=ortholog_reference),
                "mapping_confidence": conf,
                "mapping_release": mapping_release,
                "mapping_evidence": _first_present(row, ["mapping_evidence"], default=mapping_evidence),
            })

            if std_id not in seen_annotations:
                seen_annotations.add(std_id)
                a_writer.writerow({
                    "feature_id_standardized": std_id,
                    "gene_symbol": _first_present(row, ["gene_symbol"], default=orig_id),
                    "description": _first_present(row, ["description"], default=orig_id),
                    "gene_biotype": _first_present(row, ["feature_biotype", "gene_biotype"], default=feature_type),
                    "annotation_release": reference_release,
                })

    return validate_processed_table(processed_path, mapping_path, annotation_path)
