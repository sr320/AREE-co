import math

import pandas as pd

from aree.harmonize.identifiers import mapping_score
from aree.paths import root_path


# Relative importance of each scoring component. Values are normalized by their
# sum at scoring time, so the weighted evidence score always falls in [0, 1]
# regardless of the individual magnitudes chosen here.
SCORE_WEIGHTS = {
    "n_studies": 0.18,
    "total_sample_size": 0.08,
    "effect_magnitude": 0.12,
    "significance": 0.12,
    "direction_consistency": 0.14,
    "phenotype_relevance": 0.10,
    "context_breadth": 0.08,
    "assay_diversity": 0.10,
    "mapping_confidence": 0.05,
    "data_quality": 0.08,
}
_WEIGHT_TOTAL = sum(SCORE_WEIGHTS.values())


def _bounded(value, maximum):
    return min(float(value) / maximum, 1.0)


def _significance_score(q):
    if pd.isna(q):
        return 0.0
    return min(max(-math.log10(max(float(q), 1e-300)) / 10.0, 0.0), 1.0)


def score_candidates(evidence_path=None, meta_path=None, output_path=None):
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    meta_path = meta_path or root_path("data", "demo", "meta_analysis.tsv")
    evidence = pd.read_csv(evidence_path, sep="\t")
    meta = pd.read_csv(meta_path, sep="\t") if meta_path.exists() else pd.DataFrame()
    rows = []
    for feature, group in evidence.groupby("feature_id_standardized"):
        n_studies = group["study_id"].nunique()
        assay_diversity = group["feature_type"].nunique()
        total_n = group.drop_duplicates("study_id")["sample_size"].sum()
        direction_consistency = max((group["effect_size"] > 0).mean(), (group["effect_size"] < 0).mean())
        phenotype_relevance = (group["resilience_classification"] == "resilience_associated").mean()
        mapping_conf = group["mapping_confidence"].map(mapping_score).mean()
        data_quality = 1.0 - min((group["quality_flags"] != "none").mean(), 1.0) * 0.4
        best_q = group["adjusted_p_value"].min()
        effect_magnitude = min(group["effect_size"].abs().mean() / 2.5, 1.0)
        heterogeneity_penalty = 0.0
        if not meta.empty:
            match = meta[meta["feature_id_standardized"] == feature]
            if not match.empty:
                heterogeneity_penalty = min(match["i2_percent"].max() / 100.0, 1.0) * 0.15
        components = {
            "n_studies": _bounded(n_studies, 4),
            "total_sample_size": _bounded(total_n, 100),
            "effect_magnitude": effect_magnitude,
            "significance": _significance_score(best_q),
            "direction_consistency": direction_consistency,
            "phenotype_relevance": phenotype_relevance,
            "context_breadth": _bounded(group["tissue"].nunique() + group["life_stage"].nunique(), 5),
            "assay_diversity": _bounded(assay_diversity, 3),
            "mapping_confidence": mapping_conf,
            "data_quality": data_quality,
        }
        weighted = sum(SCORE_WEIGHTS[key] * components[key] for key in SCORE_WEIGHTS) / _WEIGHT_TOTAL
        score = weighted - heterogeneity_penalty
        if n_studies >= 2 and direction_consistency >= 0.67 and phenotype_relevance >= 0.5:
            category = "High-priority cross-study candidate"
        elif assay_diversity >= 2:
            category = "Multi-omics convergence candidate"
        else:
            category = "Emerging candidate requiring replication"
        if direction_consistency < 0.67:
            consistency_flag = "conflicting or context-dependent"
        else:
            consistency_flag = "reasonably consistent"
        rows.append(
            {
                "candidate_id": feature,
                "score": round(max(score, 0.0), 4),
                "category": category,
                "n_studies": n_studies,
                "total_biological_sample_size": int(total_n),
                "assay_diversity": assay_diversity,
                "direction_consistency": round(direction_consistency, 3),
                "consistency_flag": consistency_flag,
                "best_adjusted_p_value": best_q,
                "mean_mapping_confidence_score": round(mapping_conf, 3),
                "known_limitations": "; ".join(sorted(set(group["quality_flags"].astype(str)))),
            }
        )
    out = pd.DataFrame(rows).sort_values(["score", "n_studies"], ascending=False)
    output_path = output_path or root_path("data", "demo", "candidate_scores.tsv")
    out.to_csv(output_path, sep="\t", index=False)
    return output_path

