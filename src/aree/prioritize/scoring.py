import math

from pathlib import Path

import pandas as pd

from aree.groups import column_groups, is_missing, present
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


SCORED_COLUMNS = [
    "study_id",
    "feature_type",
    "sample_size",
    "effect_size",
    "effect_size_type",
    "resilience_classification",
    "mapping_confidence",
    "quality_flags",
    "adjusted_p_value",
    "tissue",
    "life_stage",
]


def _mean(values):
    return sum(values) / len(values)


def _majority_sign_fraction(effects):
    return max(sum(effect > 0 for effect in effects), sum(effect < 0 for effect in effects)) / len(effects)


def _typed_effect_summary_arrays(effects, effect_types):
    """Effect magnitude and direction consistency computed within each effect_size_type.

    Effects on different scales (e.g. log2 fold change vs methylation difference) are never
    averaged together; per-type values are combined afterwards. Direction consistency only
    counts types with at least two effects, since a lone effect cannot agree or disagree with
    anything; it is NaN when no type is replicated.
    """
    by_type = {}
    for effect, effect_type in zip(effects, effect_types):
        if not is_missing(effect_type):
            by_type.setdefault(effect_type, []).append(float(effect))
    magnitudes, directions, weights = [], [], []
    for effect_type in sorted(by_type):
        typed = by_type[effect_type]
        magnitudes.append(min(_mean([abs(effect) for effect in typed]) / 2.5, 1.0))
        if len(typed) >= 2:
            directions.append(_majority_sign_fraction(typed))
            weights.append(len(typed))
    effect_magnitude = _mean(magnitudes) if magnitudes else 0.0
    if not weights:
        return effect_magnitude, float("nan")
    direction_consistency = sum(d * w for d, w in zip(directions, weights)) / sum(weights)
    return effect_magnitude, direction_consistency


def _typed_effect_summary(group):
    return _typed_effect_summary_arrays(group["effect_size"].to_numpy(), group["effect_size_type"].to_numpy())


def _first_per_study_total(study_ids, sample_sizes):
    # Each study's sample size counts once (its first row), as drop_duplicates("study_id") did.
    seen = {}
    for study, size in zip(study_ids, sample_sizes):
        seen.setdefault(study, size)
    return sum(size for size in seen.values() if not is_missing(size))


def score_table(evidence):
    # Heterogeneity comes from the same evidence being scored, never from a previously written
    # (possibly filtered or stale) meta-analysis file. Groups with nothing poolable have no I2
    # and carry no heterogeneity information.
    meta = meta_analysis_table(evidence).dropna(subset=["i2_percent"])
    max_i2 = meta.groupby("feature_id_standardized")["i2_percent"].max().to_dict()
    rows = []
    for feature, group in column_groups(evidence, "feature_id_standardized", SCORED_COLUMNS):
        n_studies = len(set(present(group["study_id"])))
        assay_diversity = len(set(present(group["feature_type"])))
        total_n = _first_per_study_total(group["study_id"], group["sample_size"])
        effect_magnitude, direction_consistency = _typed_effect_summary_arrays(
            group["effect_size"], group["effect_size_type"]
        )
        phenotype_relevance = _mean([value == "resilience_associated" for value in group["resilience_classification"]])
        mapping_conf = _mean([mapping_score(value) for value in group["mapping_confidence"]])
        data_quality = 1.0 - min(_mean([value != "none" for value in group["quality_flags"]]), 1.0) * 0.4
        adjusted = present(group["adjusted_p_value"])
        best_q = min(adjusted) if adjusted else float("nan")
        heterogeneity_penalty = 0.0
        if feature in max_i2:
            heterogeneity_penalty = min(max_i2[feature] / 100.0, 1.0) * 0.15
        context_breadth = len(set(present(group["tissue"]))) + len(set(present(group["life_stage"])))
        components = {
            "n_studies": _bounded(n_studies, 4),
            "total_sample_size": _bounded(total_n, 100),
            "effect_magnitude": effect_magnitude,
            "significance": _significance_score(best_q),
            # Unreplicated evidence is neither rewarded nor penalized for direction.
            "direction_consistency": 0.5 if pd.isna(direction_consistency) else direction_consistency,
            "phenotype_relevance": phenotype_relevance,
            "context_breadth": _bounded(context_breadth, 5),
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
                "known_limitations": "; ".join(sorted({str(value) for value in group["quality_flags"]})),
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

