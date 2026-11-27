import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from aree.meta_analysis.random_effects import meta_analysis_table, random_effects
from aree.meta_analysis.synthesis import benjamini_hochberg, load_synthesis_config, synthesize_stratum

ROOT = Path(__file__).resolve().parents[1]


def _frame(rows):
    base = {
        "feature_id_standardized": "G1",
        "feature_type": "gene",
        "effect_size_type": "log2_fold_change",
        "phenotype": "thermal_tolerance",
        "stressor": "temperature",
        "quality_flags": "none",
    }
    return pd.DataFrame([dict(base, **row) for row in rows])


def test_block_gls_covariance_matches_explicit_matrix_solution():
    # Two shared-control contrasts (r1=2) + one 3-contrast group (r2=3) + one independent study.
    rows = [
        {"study_id": "A1", "effect_size": 1.20, "standard_error": 0.20, "quality_flags": "dependence_group=G1"},
        {"study_id": "A2", "effect_size": 1.60, "standard_error": 0.25, "quality_flags": "dependence_group=G1"},
        {"study_id": "B1", "effect_size": 0.80, "standard_error": 0.20, "quality_flags": "dependence_group=G2"},
        {"study_id": "B2", "effect_size": 1.10, "standard_error": 0.22, "quality_flags": "dependence_group=G2"},
        {"study_id": "B3", "effect_size": 0.90, "standard_error": 0.25, "quality_flags": "dependence_group=G2"},
        {"study_id": "C1", "effect_size": 0.40, "standard_error": 0.18, "quality_flags": "none"},
    ]
    rho = 0.5
    df = _frame(rows)
    res = random_effects(df, method="reml", shared_control_rho=rho)
    assert res["n_effects"] == 6
    assert res["n_studies"] == 6
    assert res["n_independent_groups"] == 3
    assert res["pooling_status"] == "pooled_covariance_adjusted"

    # Explicit full 6x6 covariance matrix Sigma = V(rho) + tau2 * Z Z^T
    sei = np.array([r["standard_error"] for r in rows])
    yi = np.array([r["effect_size"] for r in rows])
    groups = ["G1", "G1", "G2", "G2", "G2", "C1"]
    V = np.diag(sei ** 2)
    Z = np.zeros((6, 3))
    g_map = {"G1": 0, "G2": 1, "C1": 2}
    for i, gi in enumerate(groups):
        Z[i, g_map[gi]] = 1.0
        for j, gj in enumerate(groups):
            if i != j and gi == gj:
                V[i, j] = rho * sei[i] * sei[j]
    Sigma = V + res["tau2"] * (Z @ Z.T)
    ones = np.ones(6)
    inv_ones = np.linalg.solve(Sigma, ones)
    inv_y = np.linalg.solve(Sigma, yi)
    expected_mu = float((ones @ inv_y) / (ones @ inv_ones))
    expected_se = math.sqrt(1.0 / float(ones @ inv_ones))
    assert res["pooled_effect"] == pytest.approx(expected_mu, abs=1e-12)
    assert res["pooled_standard_error"] == pytest.approx(expected_se, abs=1e-12)


def test_reml_and_dl_match_metafor_reference_values():
    # Reference values verified against R metafor 4.8-0 (rma.uni & rma.mv).
    # Case 1: Independent heterogeneous k=4
    het4 = _frame([
        {"study_id": "S1", "effect_size": 0.40, "standard_error": 0.15},
        {"study_id": "S2", "effect_size": 1.25, "standard_error": 0.22},
        {"study_id": "S3", "effect_size": 0.85, "standard_error": 0.18},
        {"study_id": "S4", "effect_size": 1.60, "standard_error": 0.25},
    ])
    dl = random_effects(het4, method="dl")
    reml = random_effects(het4, method="reml")
    assert dl["q"] == pytest.approx(21.3051604861, abs=1e-8)
    assert dl["tau2"] == pytest.approx(0.2296126965, abs=1e-8)
    assert dl["pooled_effect"] == pytest.approx(1.0000582215, abs=1e-8)
    assert dl["pooled_standard_error"] == pytest.approx(0.2599079156, abs=1e-8)

    assert reml["tau2"] == pytest.approx(0.2278801788, abs=1e-8)
    assert reml["pooled_effect"] == pytest.approx(0.9998978056, abs=1e-8)
    assert reml["pooled_standard_error"] == pytest.approx(0.2590707017, abs=1e-8)
    assert reml["i2_percent"] == pytest.approx(85.8270222248, abs=1e-6)

    # Case 2: Shared-control pairs + independent study (m=3 groups, k=5 contrasts, rho=0.5)
    shared = _frame([
        {"study_id": "A1", "effect_size": 1.20, "standard_error": 0.20, "quality_flags": "dependence_group=G1"},
        {"study_id": "A2", "effect_size": 1.60, "standard_error": 0.25, "quality_flags": "dependence_group=G1"},
        {"study_id": "B1", "effect_size": -0.40, "standard_error": 0.30, "quality_flags": "dependence_group=G2"},
        {"study_id": "B2", "effect_size": -0.20, "standard_error": 0.35, "quality_flags": "dependence_group=G2"},
        {"study_id": "C1", "effect_size": 0.90, "standard_error": 0.18, "quality_flags": "none"},
    ])
    mv_reml = random_effects(shared, method="reml", shared_control_rho=0.5)
    assert mv_reml["pooled_effect"] == pytest.approx(0.6467132, abs=1e-6)
    assert mv_reml["pooled_standard_error"] == pytest.approx(0.4845994, abs=1e-6)
    assert mv_reml["tau2"] == pytest.approx(0.6567358, abs=1e-6)
    assert mv_reml["i2_percent"] == pytest.approx(93.70709, abs=1e-4)
    assert mv_reml["p_value"] == pytest.approx(0.1820297, abs=1e-6)


def test_single_independent_group_never_reports_pooled_or_replication():
    shared_only = _frame([
        {"study_id": "A1", "effect_size": 1.0, "standard_error": 0.2, "quality_flags": "dependence_group=G1"},
        {"study_id": "A2", "effect_size": 1.4, "standard_error": 0.25, "quality_flags": "dependence_group=G1"},
    ])
    # Default (rho=None) still fails closed.
    closed = random_effects(shared_only)
    assert closed["pooling_status"] == "dependent_effects_require_covariance"
    assert pd.isna(closed["pooled_effect"])

    # With rho=0.5, single independent group gives GLS within-group composite without heterogeneity/direction claims.
    comp = random_effects(shared_only, method="reml", shared_control_rho=0.5)
    assert comp["pooling_status"] == "single_independent_group"
    assert comp["n_effects"] == 2
    assert comp["n_independent_groups"] == 1
    assert pd.isna(comp["i2_percent"])
    assert pd.isna(comp["tau2"])
    assert pd.isna(comp["direction_consistency"])


def test_invalid_pooling_options_raise_value_error():
    df = _frame([{"study_id": "A", "effect_size": 1.0, "standard_error": 0.2}])
    with pytest.raises(ValueError, match="Unsupported meta-analysis method"):
        meta_analysis_table(df, method="invalid")
    for bad_rho in (-0.1, 1.0, 1.5, float("nan")):
        with pytest.raises(ValueError, match="shared_control_rho"):
            meta_analysis_table(df, shared_control_rho=bad_rho)


def test_bh_adjustment_preserves_missing_values():
    adj = benjamini_hochberg([0.01, float("nan"), 0.04, 0.02])
    assert adj[0] == pytest.approx(0.03)
    assert pd.isna(adj[1])
    assert adj[2] == pytest.approx(0.04)
    assert adj[3] == pytest.approx(0.03)


def test_synthesized_stratum_excludes_single_group_features_and_keeps_guardrails():
    ev = _frame([
        {"study_id": "A", "feature_id_standardized": "SHARED_GENE", "effect_size": 1.0, "standard_error": 0.2},
        {"study_id": "B", "feature_id_standardized": "SHARED_GENE", "effect_size": 1.5, "standard_error": 0.3},
        {"study_id": "A", "feature_id_standardized": "SINGLE_ONLY", "effect_size": 2.0, "standard_error": 0.2},
    ])
    stratum = {
        "id": "test_stratum",
        "title": "Test",
        "phenotype": "thermal_tolerance",
        "stressor": "temperature",
        "effect_size_type": "log2_fold_change",
        "study_ids": ["A", "B"],
    }
    pooled, coverage = synthesize_stratum(ev, stratum, method="reml", shared_control_rho=0.5)
    assert list(pooled["feature_id_standardized"]) == ["SHARED_GENE"]
    assert coverage["total_features_seen"] == 2
    assert coverage["multi_group_pooled_features"] == 1
    assert coverage["single_group_only_features"] == 1
    assert bool(pooled["direction_concordant"].iloc[0]) is True


def test_synthesis_strata_yaml_accounts_for_all_harmonized_studies():
    config, studies = load_synthesis_config()
    pooled_ids = [sid for s in config["strata"] for sid in s["study_ids"]]
    assert "CGIG_THERMOTOL_RNASEQ_PRJNA694496" not in pooled_ids
    assert "CGIG_OA_RNASEQ_PRJNA735889" in pooled_ids
    for s in config["strata"]:
        bioprojects = {(studies[sid].get("accessions") or {}).get("bioproject") for sid in s["study_ids"]}
        assert len(bioprojects) >= 2
