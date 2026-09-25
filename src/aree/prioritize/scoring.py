import math

from pathlib import Path

import pandas as pd

from aree.harmonize.identifiers import mapping_score
from aree.io import read_tsv
from aree.meta_analysis.random_effects import meta_analysis_table
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

SCORE_COLUMNS = [
    "candidate_id",
    "score",
    "category",
    "n_studies",
    "total_biological_sample_size",
    "assay_diversity",
    "direction_consistency",
    "consistency_flag",
    "best_adjusted_p_value",
    "mean_mapping_confidence_score",
    "known_limitations",
]


def _bounded(value, maximum):
    return min(float(value) / maximum, 1.0)


def _significance_score(q):
    if pd.isna(q):
        return 0.0
    return min(max(-math.log10(max(float(q), 1e-300)) / 10.0, 0.0), 1.0)


def _majority_sign_fraction(effects):
    return max((effects > 0).mean(), (effects < 0).mean())


def _typed_effect_summary(group):
    """Effect magnitude and direction consistency computed within each effect_size_type.

    Effects on different scales (e.g. log2 fold change vs methylation difference) are never
    averaged together; per-type values are combined afterwards. Direction consistency only
    counts types with at least two effects, since a lone effect cannot agree or disagree with
    anything; it is NaN when no type is replicated.
    """
    magnitudes, directions, weights = [], [], []
    for _, typed in group.groupby("effect_size_type"):
        magnitudes.append(min(typed["effect_size"].abs().mean() / 2.5, 1.0))
        if len(typed) >= 2:
            directions.append(_majority_sign_fraction(typed["effect_size"]))
            weights.append(len(typed))
    effect_magnitude = sum(magnitudes) / len(magnitudes)
    if not weights:
        return effect_magnitude, float("nan")
    direction_consistency = sum(d * w for d, w in zip(directions, weights)) / sum(weights)
    return effect_magnitude, direction_consistency


def score_table(evidence):
    # Heterogeneity comes from the same evidence being scored, never from a previously written
    # (possibly filtered or stale) meta-analysis file.
    meta = meta_analysis_table(evidence)
    rows = []
    for feature, group in evidence.groupby("feature_id_standardized"):
        n_studies = group["study_id"].nunique()
        assay_diversity = group["feature_type"].nunique()
        total_n = group.drop_duplicates("study_id")["sample_size"].sum()
        effect_magnitude, direction_consistency = _typed_effect_summary(group)
        phenotype_relevance = (group["resilience_classification"] == "resilience_associated").mean()
        mapping_conf = group["mapping_confidence"].map(mapping_score).mean()
        data_quality = 1.0 - min((group["quality_flags"] != "none").mean(), 1.0) * 0.4
        best_q = group["adjusted_p_value"].min()
        heterogeneity_penalty = 0.0
        # Groups with nothing poolable have no I2; they carry no heterogeneity information.
        i2 = meta.loc[meta["feature_id_standardized"] == feature, "i2_percent"].dropna()
        if not i2.empty:
            heterogeneity_penalty = min(i2.max() / 100.0, 1.0) * 0.15
        components = {
            "n_studies": _bounded(n_studies, 4),
            "total_sample_size": _bounded(total_n, 100),
            "effect_magnitude": effect_magnitude,
            "significance": _significance_score(best_q),
            # Unreplicated evidence is neither rewarded nor penalized for direction.
            "direction_consistency": 0.5 if pd.isna(direction_consistency) else direction_consistency,
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
        if pd.isna(direction_consistency):
            consistency_flag = "not replicated within an effect-size type"
        elif direction_consistency < 0.67:
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
    out = pd.DataFrame(rows, columns=SCORE_COLUMNS)
    return out.sort_values(["score", "n_studies", "candidate_id"], ascending=[False, False, True])


def score_candidates(evidence_path=None, output_path=None, phenotype=None, stressor=None):
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    evidence = read_tsv(evidence_path)
    if phenotype:
        evidence = evidence[evidence["phenotype"] == phenotype]
    if stressor:
        evidence = evidence[evidence["stressor"] == stressor]
    out = score_table(evidence)
    output_path = Path(output_path) if output_path else root_path("data", "demo", "candidate_scores.tsv")
    out.to_csv(output_path, sep="\t", index=False)
    return output_path

