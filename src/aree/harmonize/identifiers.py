import pandas as pd

from aree.io import read_tsv
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
    return read_tsv(path)


def _unresolved(original_id):
    return {
        "feature_id_standardized": original_id,
        "ortholog_reference": None,
        "mapping_confidence": "unresolved",
    }


def _mapping_result(hit):
    ortholog_reference = hit.get("ortholog_reference")
    if pd.isna(ortholog_reference):
        ortholog_reference = None
    result = {
        "feature_id_standardized": hit["feature_id_standardized"],
        "ortholog_reference": ortholog_reference,
        "mapping_confidence": hit["mapping_confidence"],
    }
    for field in ("mapping_release", "mapping_evidence"):
        value = hit.get(field)
        if value is not None and not pd.isna(value):
            result[field] = value
    return result


def mapping_lookup(mapping):
    """Index a mapping table by original ID; the first row wins, as in ``map_identifier``."""
    lookup = {}
    for hit in mapping.to_dict("records"):
        original_id = hit["feature_id_original"]
        if not pd.isna(original_id):
            lookup.setdefault(original_id, hit)
    return lookup


def map_identifiers(original_ids, mapping=None):
    """Map many IDs with one pass over the mapping table (``map_identifier`` scans it per ID)."""
    lookup = mapping_lookup(mapping if mapping is not None else load_mapping())
    return [
        _mapping_result(lookup[original_id]) if original_id in lookup else _unresolved(original_id)
        for original_id in original_ids
    ]


def map_identifier(original_id, mapping=None):
    mapping = mapping if mapping is not None else load_mapping()
    hits = mapping[mapping["feature_id_original"] == original_id]
    if hits.empty:
        return _unresolved(original_id)
    return _mapping_result(hits.iloc[0].to_dict())


def mapping_score(confidence):
    return CONFIDENCE_ORDER.get(confidence, 0) / 5.0
