import hashlib
from datetime import date
from pathlib import Path

import pandas as pd

from aree import __version__
from aree.harmonize.identifiers import load_mapping, map_identifier
from aree.intake.registry import load_study
from aree.paths import root_path
from aree.validation.schemas import validate_evidence_records


REQUIRED_PROCESSED_COLUMNS = [
    "sample_comparison",
    "feature_id_original",
    "feature_type",
    "molecular_direction",
    "effect_size",
    "effect_size_type",
    "standard_error",
    "p_value",
    "adjusted_p_value",
    "analysis_method",
]


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def harmonize_processed(study_id, input_path, output_path=None):
    study = load_study(study_id)
    table = pd.read_csv(input_path, sep="\t")
    missing = [column for column in REQUIRED_PROCESSED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError("Processed result table is missing columns: {}".format(", ".join(missing)))

    mapping = load_mapping()
    input_checksum = checksum(input_path)
    records = []
    for idx, row in table.iterrows():
        mapped = map_identifier(row["feature_id_original"], mapping)
        record = {
            "evidence_id": "{}:{}".format(study_id, idx + 1),
            "study_id": study_id,
            "sample_comparison": row["sample_comparison"],
            "feature_id_original": str(row["feature_id_original"]),
            "feature_id_standardized": str(mapped["feature_id_standardized"]),
            "feature_type": row["feature_type"],
            "species": study["species"],
            "genome_assembly": study.get("genome_assembly"),
            "annotation_version": study.get("annotation_version"),
            "ortholog_reference": mapped["ortholog_reference"],
            "mapping_confidence": mapped["mapping_confidence"],
            "molecular_direction": row["molecular_direction"],
            "effect_size": float(row["effect_size"]),
            "effect_size_type": row["effect_size_type"],
            "standard_error": None if pd.isna(row["standard_error"]) else float(row["standard_error"]),
            "ci_lower": None,
            "ci_upper": None,
            "p_value": None if pd.isna(row["p_value"]) else float(row["p_value"]),
            "adjusted_p_value": None if pd.isna(row["adjusted_p_value"]) else float(row["adjusted_p_value"]),
            "sample_size": int(study["sample_size"]),
            "tissue": study["tissue"],
            "life_stage": study["life_stage"],
            "stressor": study["stressor_class"],
            "phenotype": study["phenotype"],
            "phenotype_direction": study["phenotype_direction"],
            "analysis_method": row["analysis_method"],
            "quality_flags": row.get("quality_flags", "none"),
            "source_file": str(input_path),
            "input_checksum": input_checksum,
            "workflow_version": __version__,
            "date_generated": date.today().isoformat(),
            "source_accession": ";".join(
                str(v) for v in study.get("accessions", {}).values() if v is not None
            ),
            "resilience_classification": study["resilience_classification"],
        }
        records.append(record)
    validate_evidence_records(records)
    output_path = Path(output_path) if output_path else root_path("data", "demo", "harmonized_evidence.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    new_table = pd.DataFrame(records)
    if output_path.exists() and output_path.stat().st_size > 0:
        existing = pd.read_csv(output_path, sep="\t")
        existing = existing[existing["study_id"] != study_id]
        new_table = pd.concat([existing, new_table], ignore_index=True)
    new_table.to_csv(output_path, sep="\t", index=False)
    return output_path


def harmonize_demo():
    processed_dir = root_path("data", "demo", "processed")
    for path in sorted(processed_dir.glob("*_*.tsv")):
        study_id = path.name.rsplit("_", 1)[0]
        harmonize_processed(study_id, path)
    return root_path("data", "demo", "harmonized_evidence.tsv")

