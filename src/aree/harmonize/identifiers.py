import pandas as pd

from aree.paths import root_path


CONFIDENCE_ORDER = {
    "exact": 5,
    "one-to-one ortholog": 4,
    "many-to-one ortholog": 3,
    "one-to-many ortholog": 2,
    "inferred": 1,
    "unresolved": 0,
}


def load_mapping(path=None):
    path = path or root_path("data", "mappings", "demo_identifier_map.tsv")
    return pd.read_csv(path, sep="\t")


def map_identifier(original_id, mapping=None):
    mapping = mapping if mapping is not None else load_mapping()
    hits = mapping[mapping["feature_id_original"] == original_id]
    if hits.empty:
        return {
            "feature_id_standardized": original_id,
            "ortholog_reference": None,
            "mapping_confidence": "unresolved",
        }
    hit = hits.iloc[0].to_dict()
    ortholog_reference = hit.get("ortholog_reference")
    if pd.isna(ortholog_reference):
        ortholog_reference = None
    return {
        "feature_id_standardized": hit["feature_id_standardized"],
        "ortholog_reference": ortholog_reference,
        "mapping_confidence": hit["mapping_confidence"],
    }


def mapping_score(confidence):
    return CONFIDENCE_ORDER.get(confidence, 0) / 5.0
