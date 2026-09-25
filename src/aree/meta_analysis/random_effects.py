import math
from pathlib import Path

import pandas as pd

from aree.groups import column_groups, is_missing
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


def _pool(effects, standard_errors, study_ids):
    """DerSimonian-Laird pooling of one group's arrays, reporting effects that could not be pooled.

    Sums are plain sequential sums, which are deterministic on every platform.
    """
    usable = [
        not is_missing(effect) and not is_missing(se) and se > 0
        for effect, se in zip(effects, standard_errors)
    ]
    excluded = sorted({study for study, ok in zip(study_ids, usable) if not ok})
    exclusion = {
        "n_effects_excluded": usable.count(False),
        "excluded_study_ids": ";".join(excluded),
    }
    yi = [float(effect) for effect, ok in zip(effects, usable) if ok]
    # x * x, not x ** 2: pow() is not exactly rounded on every platform, multiplication is.
    vi = [float(se) * float(se) for se, ok in zip(standard_errors, usable) if ok]
    studies = sorted({study for study, ok in zip(study_ids, usable) if ok})
    k = len(yi)
    if k == 0:
        # Keep the group visible: silently dropping it hides studies that lack standard errors.
        result = {column: float("nan") for column in RESULT_COLUMNS}
        result.update(n_effects=0, n_studies=0, study_ids="", pooling_status="no_standard_errors", **exclusion)
        return result
    wi = [1.0 / v for v in vi]
    sum_w = sum(wi)
    fixed = sum(w * y for w, y in zip(wi, yi)) / sum_w
    q = sum(w * (y - fixed) * (y - fixed) for w, y in zip(wi, yi))
    c = sum_w - sum(w * w for w in wi) / sum_w
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    rei = [1.0 / (v + tau2) for v in vi]
    sum_re = sum(rei)
    pooled = sum(r * y for r, y in zip(rei, yi)) / sum_re
    se = math.sqrt(1.0 / sum_re)
    z = pooled / se if se > 0 else 0.0
    i2 = max(0.0, (q - (k - 1)) / q) * 100.0 if q > 0 and k > 1 else 0.0
    return {
        "n_effects": k,
        "n_studies": len(studies),
        "pooled_effect": pooled,
        "pooled_standard_error": se,
        "p_value": _two_sided_p(z),
        "q": q,
        "i2_percent": i2,
        "tau2": tau2,
        "direction_consistency": max(sum(y > 0 for y in yi), sum(y < 0 for y in yi)) / k,
        "study_ids": ";".join(studies),
        "pooling_status": "pooled" if k > 1 else "single_effect",
        **exclusion,
    }


def random_effects(group):
    """DerSimonian-Laird pooling of one group, reporting any effects that could not be pooled."""
    return _pool(
        group["effect_size"].to_numpy(), group["standard_error"].to_numpy(), group["study_id"].to_numpy()
    )


GROUP_COLUMNS = ["feature_id_standardized", "feature_type", "effect_size_type", "phenotype", "stressor"]


def meta_analysis_table(evidence):
    """Pool effects per feature and context; only effects on the same scale (effect_size_type) are pooled."""
    rows = []
    groups = column_groups(evidence, GROUP_COLUMNS, ["effect_size", "standard_error", "study_id"])
    for keys, columns in groups:
        row = dict(zip(GROUP_COLUMNS, keys))
        row.update(_pool(columns["effect_size"], columns["standard_error"], columns["study_id"]))
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
