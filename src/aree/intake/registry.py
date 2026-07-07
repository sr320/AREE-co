from pathlib import Path

import pandas as pd

from aree.paths import root_path
from aree.validation.schemas import validate_study_file


REGISTRY_COLUMNS = [
    "study_id",
    "doi",
    "citation",
    "species",
    "assay_type",
    "tissue",
    "life_stage",
    "stressor_class",
    "phenotype",
    "resilience_classification",
    "sample_size",
    "data_status",
    "analysis_status",
    "quality_control_status",
]


def flatten_for_registry(study):
    return {
        "study_id": study["study_id"],
        "doi": study.get("doi"),
        "citation": study["citation"],
        "species": study["species"],
        "assay_type": study["assay_type"],
        "tissue": study["tissue"],
        "life_stage": study["life_stage"],
        "stressor_class": study["stressor_class"],
        "phenotype": study["phenotype"],
        "resilience_classification": study["resilience_classification"],
        "sample_size": study["sample_size"],
        "data_status": study["data_availability"]["status"],
        "analysis_status": study["analysis_status"],
        "quality_control_status": study["quality_control_status"],
    }


def register_study(path, registry_path=None):
    study = validate_study_file(path)
    registry_path = Path(registry_path) if registry_path else root_path("registry", "study_registry.csv")
    if registry_path.exists() and registry_path.stat().st_size > 0:
        registry = pd.read_csv(registry_path)
    else:
        registry = pd.DataFrame(columns=REGISTRY_COLUMNS)
    if "study_id" in registry.columns and study["study_id"] in set(registry["study_id"].dropna()):
        raise ValueError("Duplicate study_id: {}".format(study["study_id"]))
    registry = pd.concat([registry, pd.DataFrame([flatten_for_registry(study)])], ignore_index=True)
    registry = registry[REGISTRY_COLUMNS]
    registry.to_csv(registry_path, index=False)
    return registry_path


def load_study(study_id):
    path = root_path("registry", "studies", "{}.yaml".format(study_id))
    return validate_study_file(path)

