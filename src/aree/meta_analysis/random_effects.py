import math
from pathlib import Path

import pandas as pd

from aree.groups import column_groups, dependence_groups, is_missing
from aree.io import read_evidence
from aree.paths import root_path


RESULT_COLUMNS = [
    "n_effects",
    "n_studies",
    "n_independent_groups",
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
VALID_METHODS = ("dl", "reml")


def _two_sided_p(z):
    # erfc keeps precision for large |z| (1 - cdf underflows to 0). Rounding to 12 significant
    # digits absorbs last-ulp differences between platform libm implementations, so written
    # outputs are identical on macOS and Linux.
    return float("{:.12g}".format(math.erfc(abs(z) / math.sqrt(2.0))))


EXPLORATORY_FLAG = "exploratory_tank_clustering_unmodeled"


def inference_eligible(evidence):
    if "quality_flags" not in evidence:
        return pd.Series(True, index=evidence.index)
    return ~evidence["quality_flags"].fillna("").str.split(";").map(
        lambda flags: EXPLORATORY_FLAG in flags
    ).astype(bool)


def poolable(evidence):
    """Effects with usable uncertainty and eligible study-level inference."""
    return (evidence["effect_size"].notna() & evidence["standard_error"].notna()
            & (evidence["standard_error"] > 0) & inference_eligible(evidence))


def _validate_pooling_options(method, shared_control_rho):
    if method not in VALID_METHODS:
        raise ValueError("Unsupported meta-analysis method {!r}; expected one of {}".format(method, VALID_METHODS))
    if shared_control_rho is not None:
        rho = float(shared_control_rho)
        if not (0.0 <= rho < 1.0) or math.isnan(rho):
            raise ValueError("shared_control_rho must be in [0.0, 1.0), got {!r}".format(shared_control_rho))
        return rho
    return None


def _solve_linear_system(matrix, rhs):
    """Solve ``matrix @ x = rhs`` for a small SPD covariance block via Gaussian elimination."""
    n = len(rhs)
    aug = [list(map(float, matrix[i])) + [float(rhs[i])] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) <= 0.0:
            raise ValueError("Singular sampling covariance matrix")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        diag = aug[col][col]
        for row in range(col + 1, n):
            factor = aug[row][col] / diag
            if factor != 0.0:
                for k in range(col, n + 1):
                    aug[row][k] -= factor * aug[col][k]
    sol = [0.0] * n
    for i in range(n - 1, -1, -1):
        rem = aug[i][n] - sum(aug[i][j] * sol[j] for j in range(i + 1, n))
        sol[i] = rem / aug[i][i]
    return sol


def _collapse_dependent_groups(effects, standard_errors, groups, rho):
    """Reduce each dependence group to its exact GLS summary ``(y_bar_g, v_tilde_g)``.

    For a group with ``r`` contrasts sharing controls with sampling correlation ``rho``,
    the ``r x r`` sampling covariance block ``V_g`` has ``V_aa = SE_a^2`` and
    ``V_ab = rho * SE_a * SE_b``. Under ``Sigma(tau^2) = V + tau^2 Z Z^T``, GLS
    and REML over all ``k`` contrasts are mathematically equivalent to univariate
    random-effects synthesis over these ``m`` independent group summaries.
    """
    grouped = {}
    for effect, se, group in zip(effects, standard_errors, groups):
        grouped.setdefault(group, []).append((effect, se))
    summary_y, summary_v, group_sizes = [], [], []
    for group in sorted(grouped):
        items = grouped[group]
        r = len(items)
        group_sizes.append(r)
        if r == 1:
            y0, se0 = items[0]
            summary_y.append(y0)
            summary_v.append(se0 * se0)
        else:
            cov = [
                [
                    (items[i][1] * items[i][1]) if i == j else (rho * items[i][1] * items[j][1])
                    for j in range(r)
                ]
                for i in range(r)
            ]
            ones = [1.0] * r
            ys = [item[0] for item in items]
            v_inv_1 = _solve_linear_system(cov, ones)
            v_inv_y = _solve_linear_system(cov, ys)
            precision = sum(v_inv_1)
            summary_y.append(sum(v_inv_y) / precision)
            summary_v.append(1.0 / precision)
    return summary_y, summary_v, group_sizes


def _reml_loglik(yi, vi, tau2):
    wi = [1.0 / (v + tau2) for v in vi]
    s0 = sum(wi)
    mu = sum(w * y for w, y in zip(wi, yi)) / s0
    rss = sum(w * (y - mu) * (y - mu) for w, y in zip(wi, yi))
    logdet = sum(math.log(v + tau2) for v in vi) + math.log(s0)
    return -0.5 * (logdet + rss)


def _reml_tau2(yi, vi, tau2_init):
    """Exact REML heterogeneity variance via safeguarded Fisher scoring."""
    tau2 = max(0.0, float(tau2_init))
    ll_cur = _reml_loglik(yi, vi, tau2)
    for _ in range(200):
        wi = [1.0 / (v + tau2) for v in vi]
        s0 = sum(wi)
        mu = sum(w * y for w, y in zip(wi, yi)) / s0
        w2 = [w * w for w in wi]
        sum_w2 = sum(w2)
        sum_w3 = sum(w * w2_i for w, w2_i in zip(wi, w2))
        tr_p = s0 - sum_w2 / s0
        tr_p2 = sum_w2 - 2.0 * sum_w3 / s0 + (sum_w2 * sum_w2) / (s0 * s0)
        yp2y = sum(w2_i * (y - mu) * (y - mu) for w2_i, y in zip(w2, yi))
        score = 0.5 * (yp2y - tr_p)
        info = 0.5 * tr_p2
        if tau2 == 0.0 and score <= 0.0:
            break
        if info <= 0.0:
            break
        step = score / info
        if abs(step) <= 1e-12 * (1.0 + tau2):
            break
        alpha = 1.0
        accepted = False
        while alpha > 1e-6:
            cand = max(0.0, tau2 + alpha * step)
            ll_cand = _reml_loglik(yi, vi, cand)
            if ll_cand >= ll_cur - 1e-14:
                tau2 = cand
                ll_cur = ll_cand
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break
    if _reml_loglik(yi, vi, 0.0) > ll_cur:
        return 0.0
    return tau2


def _pool(effects, standard_errors, study_ids, eligible=None, replication_groups=None, method="dl", shared_control_rho=None):
    """Random-effects pooling (DL or REML), with optional shared-control covariance modeling.

    When ``shared_control_rho`` is ``None`` and multiple usable effects share a dependence group,
    pooling fails closed with ``dependent_effects_require_covariance``. When ``shared_control_rho``
    is provided in ``[0, 1)``, contrasts sharing a dependence group are collapsed via exact
    block-covariance GLS before random-effects synthesis across independent groups.
    """
    rho = _validate_pooling_options(method, shared_control_rho)
    usable = [
        not is_missing(effect) and not is_missing(se) and se > 0
        for effect, se in zip(effects, standard_errors)
    ]
    if eligible is not None:
        usable = [ok and bool(allowed) for ok, allowed in zip(usable, eligible)]
    if replication_groups is not None:
        included_groups = [group for group, ok in zip(replication_groups, usable) if ok]
        if rho is None and len(set(included_groups)) < len(included_groups):
            # Shared samples require covariances that the evidence table does not supply.
            # Retain the whole context for review rather than choose an arbitrary contrast.
            result = {column: float("nan") for column in RESULT_COLUMNS}
            result.update(n_effects=0, n_studies=0, n_independent_groups=0, study_ids="",
                          pooling_status="dependent_effects_require_covariance",
                          n_effects_excluded=len(effects), excluded_study_ids=";".join(sorted(set(study_ids))))
            return result
    excluded = sorted({study for study, ok in zip(study_ids, usable) if not ok})
    exclusion = {
        "n_effects_excluded": usable.count(False),
        "excluded_study_ids": ";".join(excluded),
    }
    yi = [float(effect) for effect, ok in zip(effects, usable) if ok]
    sei = [float(se) for se, ok in zip(standard_errors, usable) if ok]
    # x * x, not x ** 2: pow() is not exactly rounded on every platform, multiplication is.
    vi = [se * se for se in sei]
    studies = sorted({study for study, ok in zip(study_ids, usable) if ok})
    k = len(yi)
    if k == 0:
        # Keep the group visible: silently dropping it hides studies that lack standard errors.
        result = {column: float("nan") for column in RESULT_COLUMNS}
        status = "no_inference_eligible_effects" if eligible is not None and not all(eligible) else "no_standard_errors"
        result.update(n_effects=0, n_studies=0, n_independent_groups=0, study_ids="", pooling_status=status, **exclusion)
        return result

    if replication_groups is not None and rho is not None:
        usable_groups = [group for group, ok in zip(replication_groups, usable) if ok]
        group_y, group_v, group_sizes = _collapse_dependent_groups(yi, sei, usable_groups, rho)
    else:
        usable_groups = studies if replication_groups is None else [
            group for group, ok in zip(replication_groups, usable) if ok
        ]
        group_y, group_v, group_sizes = yi, vi, [1] * k

    m = len(group_y)
    has_shared_group = any(size > 1 for size in group_sizes)
    if m == 1:
        pooled = group_y[0]
        se = math.sqrt(group_v[0])
        z = pooled / se if se > 0 else 0.0
        if has_shared_group:
            # Correlated contrasts within a single replication group cannot establish
            # between-study heterogeneity or independent directional replication.
            return {
                "n_effects": k,
                "n_studies": len(studies),
                "n_independent_groups": 1,
                "pooled_effect": pooled,
                "pooled_standard_error": se,
                "p_value": _two_sided_p(z),
                "q": float("nan"),
                "i2_percent": float("nan"),
                "tau2": float("nan"),
                "direction_consistency": float("nan"),
                "study_ids": ";".join(studies),
                "pooling_status": "single_independent_group",
                **exclusion,
            }
        return {
            "n_effects": 1,
            "n_studies": len(studies),
            "n_independent_groups": 1,
            "pooled_effect": pooled,
            "pooled_standard_error": se,
            "p_value": _two_sided_p(z),
            "q": 0.0,
            "i2_percent": 0.0,
            "tau2": 0.0,
            "direction_consistency": 1.0,
            "study_ids": ";".join(studies),
            "pooling_status": "single_effect",
            **exclusion,
        }

    wi = [1.0 / v for v in group_v]
    sum_w = sum(wi)
    fixed = sum(w * y for w, y in zip(wi, group_y)) / sum_w
    q = sum(w * (y - fixed) * (y - fixed) for w, y in zip(wi, group_y))
    c = sum_w - sum(w * w for w in wi) / sum_w
    tau2_dl = max(0.0, (q - (m - 1)) / c) if c > 0 else 0.0
    if method == "reml":
        tau2 = _reml_tau2(group_y, group_v, tau2_dl)
        v_typ = (m - 1) / c if c > 0 else 0.0
        i2 = (tau2 / (tau2 + v_typ)) * 100.0 if tau2 > 0 and v_typ > 0 else 0.0
    else:
        tau2 = tau2_dl
        i2 = max(0.0, (q - (m - 1)) / q) * 100.0 if q > 0 else 0.0
    rei = [1.0 / (v + tau2) for v in group_v]
    sum_re = sum(rei)
    pooled = sum(r * y for r, y in zip(rei, group_y)) / sum_re
    se = math.sqrt(1.0 / sum_re)
    z = pooled / se if se > 0 else 0.0
    return {
        "n_effects": k,
        "n_studies": len(studies),
        "n_independent_groups": m,
        "pooled_effect": pooled,
        "pooled_standard_error": se,
        "p_value": _two_sided_p(z),
        "q": q,
        "i2_percent": i2,
        "tau2": tau2,
        "direction_consistency": max(sum(y > 0 for y in group_y), sum(y < 0 for y in group_y)) / m,
        "study_ids": ";".join(studies),
        "pooling_status": "pooled_covariance_adjusted" if has_shared_group else "pooled",
        **exclusion,
    }


def random_effects(group, method="dl", shared_control_rho=None):
    """Random-effects pooling of one group, reporting any effects that could not be pooled."""
    return _pool(
        group["effect_size"].to_numpy(), group["standard_error"].to_numpy(), group["study_id"].to_numpy(),
        inference_eligible(group).to_numpy(),
        dependence_groups(group),
        method=method,
        shared_control_rho=shared_control_rho,
    )


GROUP_COLUMNS = ["feature_id_standardized", "feature_type", "effect_size_type", "phenotype", "stressor"]


def meta_analysis_table(evidence, method="dl", shared_control_rho=None):
    """Pool effects per feature and context; only effects on the same scale (effect_size_type) are pooled."""
    _validate_pooling_options(method, shared_control_rho)
    rows = []
    evidence = evidence.assign(_inference_eligible=inference_eligible(evidence),
                               _replication_group=dependence_groups(evidence))
    groups = column_groups(evidence, GROUP_COLUMNS, ["effect_size", "standard_error", "study_id", "_inference_eligible", "_replication_group"])
    for keys, columns in groups:
        row = dict(zip(GROUP_COLUMNS, keys))
        row.update(_pool(
            columns["effect_size"], columns["standard_error"], columns["study_id"],
            columns["_inference_eligible"], columns["_replication_group"],
            method=method, shared_control_rho=shared_control_rho,
        ))
        rows.append(row)
    return pd.DataFrame(rows, columns=GROUP_COLUMNS + RESULT_COLUMNS)


def run_meta_analysis(
    phenotype=None, feature_type=None, evidence_path=None, output_path=None,
    method="dl", shared_control_rho=None,
):
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    evidence = read_evidence(evidence_path)
    if phenotype:
        evidence = evidence[evidence["phenotype"] == phenotype]
    if feature_type:
        evidence = evidence[evidence["feature_type"] == feature_type]
    output_path = Path(output_path) if output_path else root_path("data", "demo", "meta_analysis.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    meta_analysis_table(evidence, method=method, shared_control_rho=shared_control_rho).to_csv(output_path, sep="\t", index=False)
    return output_path
