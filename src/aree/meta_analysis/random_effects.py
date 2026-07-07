import math

import pandas as pd

from aree.paths import root_path


def _normal_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def random_effects(group):
    group = group.dropna(subset=["effect_size", "standard_error"]).copy()
    group = group[group["standard_error"] > 0]
    k = len(group)
    if k == 0:
        return None
    yi = group["effect_size"].astype(float)
    vi = group["standard_error"].astype(float) ** 2
    wi = 1.0 / vi
    fixed = (wi * yi).sum() / wi.sum()
    q = (wi * (yi - fixed) ** 2).sum()
    c = wi.sum() - (wi**2).sum() / wi.sum()
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    rei = 1.0 / (vi + tau2)
    pooled = (rei * yi).sum() / rei.sum()
    se = math.sqrt(1.0 / rei.sum())
    z = pooled / se if se > 0 else 0.0
    p_value = 2.0 * (1.0 - _normal_cdf(abs(z)))
    i2 = max(0.0, (q - (k - 1)) / q) * 100.0 if q > 0 and k > 1 else 0.0
    return {
        "n_effects": k,
        "n_studies": group["study_id"].nunique(),
        "pooled_effect": pooled,
        "pooled_standard_error": se,
        "p_value": p_value,
        "q": q,
        "i2_percent": i2,
        "tau2": tau2,
        "direction_consistency": max((yi > 0).mean(), (yi < 0).mean()),
        "study_ids": ";".join(sorted(group["study_id"].unique())),
    }


def run_meta_analysis(phenotype=None, feature_type=None, evidence_path=None, output_path=None):
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    evidence = pd.read_csv(evidence_path, sep="\t")
    if phenotype:
        evidence = evidence[evidence["phenotype"] == phenotype]
    if feature_type:
        evidence = evidence[evidence["feature_type"] == feature_type]
    rows = []
    group_cols = ["feature_id_standardized", "feature_type", "phenotype", "stressor"]
    for keys, group in evidence.groupby(group_cols):
        result = random_effects(group)
        if result:
            row = dict(zip(group_cols, keys))
            row.update(result)
            rows.append(row)
    out = pd.DataFrame(rows)
    output_path = output_path or root_path("data", "demo", "meta_analysis.tsv")
    out.to_csv(output_path, sep="\t", index=False)
    return output_path

