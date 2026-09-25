import math
from pathlib import Path

import pandas as pd

from aree.io import read_tsv
from aree.paths import root_path


RESULT_COLUMNS = [
    "n_effects",
    "n_studies",
    "pooled_effect",
    "pooled_standard_error",
    "p_value",
    "q",
    "i2_percent",
    "tau2",
    "direction_consistency",
    "study_ids",
    "pooling_status",
    "n_effects_excluded",
    "excluded_study_ids",
]


def _two_sided_p(z):
    # erfc keeps precision for large |z| (1 - cdf underflows to 0). Rounding to 12 significant
    # digits absorbs last-ulp differences between platform libm implementations, so written
    # outputs are identical on macOS and Linux.
    return float("{:.12g}".format(math.erfc(abs(z) / math.sqrt(2.0))))


def poolable(evidence):
    """Mask of effects usable for inverse-variance pooling (effect size and a positive standard error)."""
    return evidence["effect_size"].notna() & evidence["standard_error"].notna() & (evidence["standard_error"] > 0)


def random_effects(group):
    """DerSimonian-Laird pooling of one group, reporting any effects that could not be pooled."""
    usable = poolable(group)
    excluded = group[~usable]
    exclusion = {
        "n_effects_excluded": len(excluded),
        "excluded_study_ids": ";".join(sorted(excluded["study_id"].unique())),
    }
    group = group[usable]
    k = len(group)
    if k == 0:
        # Keep the group visible: silently dropping it hides studies that lack standard errors.
        result = {column: float("nan") for column in RESULT_COLUMNS}
        result.update(n_effects=0, n_studies=0, study_ids="", pooling_status="no_standard_errors", **exclusion)
        return result
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
    p_value = _two_sided_p(z)
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
        "pooling_status": "pooled" if k > 1 else "single_effect",
        **exclusion,
    }


GROUP_COLUMNS = ["feature_id_standardized", "feature_type", "effect_size_type", "phenotype", "stressor"]


def meta_analysis_table(evidence):
    """Pool effects per feature and context; only effects on the same scale (effect_size_type) are pooled."""
    rows = []
    for keys, group in evidence.groupby(GROUP_COLUMNS):
        row = dict(zip(GROUP_COLUMNS, keys))
        row.update(random_effects(group))
        rows.append(row)
    return pd.DataFrame(rows, columns=GROUP_COLUMNS + RESULT_COLUMNS)


def run_meta_analysis(phenotype=None, feature_type=None, evidence_path=None, output_path=None):
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    evidence = read_tsv(evidence_path)
    if phenotype:
        evidence = evidence[evidence["phenotype"] == phenotype]
    if feature_type:
        evidence = evidence[evidence["feature_type"] == feature_type]
    output_path = Path(output_path) if output_path else root_path("data", "demo", "meta_analysis.tsv")
    meta_analysis_table(evidence).to_csv(output_path, sep="\t", index=False)
    return output_path
