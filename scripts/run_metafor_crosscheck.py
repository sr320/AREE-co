#!/usr/bin/env python3
"""Cross-check Python DL and REML covariance-aware pooling against R metafor."""

import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from aree.io import read_evidence
from aree.meta_analysis.random_effects import random_effects

ROOT = Path(__file__).resolve().parents[1]


def _build_synthetic_cases():
    cases = []
    # 1. Homogeneous (tau2 = 0)
    cases.append(("syn_homog_dl", "dl", None, pd.DataFrame([
        {"study_id": "S1", "effect_size": 1.00, "standard_error": 0.40},
        {"study_id": "S2", "effect_size": 1.15, "standard_error": 0.45},
        {"study_id": "S3", "effect_size": 0.90, "standard_error": 0.35},
    ])))
    cases.append(("syn_homog_reml", "reml", None, cases[-1][3]))
    # 2. Moderate heterogeneity k=4
    het4 = pd.DataFrame([
        {"study_id": "S1", "effect_size": 0.40, "standard_error": 0.15},
        {"study_id": "S2", "effect_size": 1.25, "standard_error": 0.22},
        {"study_id": "S3", "effect_size": 0.85, "standard_error": 0.18},
        {"study_id": "S4", "effect_size": 1.60, "standard_error": 0.25},
    ])
    cases.append(("syn_het4_dl", "dl", None, het4))
    cases.append(("syn_het4_reml", "reml", None, het4))
    # 3. Opposing directions k=3
    opp3 = pd.DataFrame([
        {"study_id": "S1", "effect_size": -1.10, "standard_error": 0.25},
        {"study_id": "S2", "effect_size": 0.95, "standard_error": 0.30},
        {"study_id": "S3", "effect_size": 0.20, "standard_error": 0.20},
    ])
    cases.append(("syn_opp3_dl", "dl", None, opp3))
    cases.append(("syn_opp3_reml", "reml", None, opp3))
    # 4. Large z precision check
    largez = pd.DataFrame([
        {"study_id": "S1", "effect_size": 5.5, "standard_error": 0.35},
        {"study_id": "S2", "effect_size": 5.8, "standard_error": 0.40},
        {"study_id": "S3", "effect_size": 5.2, "standard_error": 0.30},
    ])
    cases.append(("syn_largez_reml", "reml", None, largez))
    # 5. Shared-control block covariance (rho=0.5 and rho=0.8) + independent study
    shared = pd.DataFrame([
        {"study_id": "A1", "effect_size": 1.20, "standard_error": 0.20, "quality_flags": "dependence_group=G1"},
        {"study_id": "A2", "effect_size": 1.60, "standard_error": 0.25, "quality_flags": "dependence_group=G1"},
        {"study_id": "B1", "effect_size": -0.40, "standard_error": 0.30, "quality_flags": "dependence_group=G2"},
        {"study_id": "B2", "effect_size": -0.20, "standard_error": 0.35, "quality_flags": "dependence_group=G2"},
        {"study_id": "C1", "effect_size": 0.90, "standard_error": 0.18, "quality_flags": "none"},
    ])
    cases.append(("syn_shared_rho05_dl", "dl", 0.5, shared))
    cases.append(("syn_shared_rho05_reml", "reml", 0.5, shared))
    cases.append(("syn_shared_rho08_reml", "reml", 0.8, shared))
    # 6. 3-contrast shared group + 1 independent
    shared3 = pd.DataFrame([
        {"study_id": "A1", "effect_size": 0.8, "standard_error": 0.20, "quality_flags": "dependence_group=G1"},
        {"study_id": "A2", "effect_size": 1.1, "standard_error": 0.22, "quality_flags": "dependence_group=G1"},
        {"study_id": "A3", "effect_size": 0.9, "standard_error": 0.25, "quality_flags": "dependence_group=G1"},
        {"study_id": "B1", "effect_size": 1.7, "standard_error": 0.19, "quality_flags": "none"},
    ])
    cases.append(("syn_shared3_rho05_reml", "reml", 0.5, shared3))
    cases.append(("syn_shared3_rho05_dl", "dl", 0.5, shared3))
    return cases


def _build_real_cases():
    ev = read_evidence(ROOT / "data" / "harmonized")
    cases = []
    # Independent-only stratum: thermal protection (PRJEB86646 + PRJNA593309)
    sub_ind = ev[ev["study_id"].isin(["CGIG_HEAT_RNASEQ_PRJEB86646", "CGIG_HEAT_RNASEQ_PRJNA593309"])]
    genes_ind = sorted(sub_ind.groupby("feature_id_standardized").filter(lambda g: len(g) == 2)["feature_id_standardized"].unique())
    step_ind = max(1, len(genes_ind) // 40)
    for gene in genes_ind[::step_ind][:40]:
        grp = sub_ind[sub_ind["feature_id_standardized"] == gene]
        cases.append(("real_ind_dl_" + gene, "dl", None, grp))
        cases.append(("real_ind_reml_" + gene, "reml", None, grp))

    # Shared-control stratum: acute heat gill (PRJNA407831 6h+24h, PRJNA516762, PRJNA913164)
    heat_studies = [
        "CGIG_HEAT6H_RNASEQ_PRJNA407831",
        "CGIG_HEAT24H_RNASEQ_PRJNA407831",
        "CGIG_HEAT_RNASEQ_PRJNA516762",
        "CGIG_HEAT_RNASEQ_PRJNA913164",
    ]
    sub_heat = ev[ev["study_id"].isin(heat_studies)]
    genes_heat = sorted(sub_heat.groupby("feature_id_standardized").filter(lambda g: len(g) == 4)["feature_id_standardized"].unique())
    step_heat = max(1, len(genes_heat) // 30)
    for gene in genes_heat[::step_heat][:30]:
        grp = sub_heat[sub_heat["feature_id_standardized"] == gene]
        cases.append(("real_heat_rho05_reml_" + gene, "reml", 0.5, grp))
        cases.append(("real_heat_rho05_dl_" + gene, "dl", 0.5, grp))
    return cases


def run_crosscheck(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = _build_synthetic_cases() + _build_real_cases()

    r_input_rows = []
    py_results = {}
    for cid, method, rho, df in cases:
        res = random_effects(df, method=method, shared_control_rho=rho)
        py_results[cid] = (method, rho, res)
        from aree.groups import dependence_groups
        dep = dependence_groups(df)
        for (_, row), g in zip(df.iterrows(), dep):
            r_input_rows.append({
                "case_id": cid,
                "method": method,
                "rho": "" if rho is None else rho,
                "study_id": row["study_id"],
                "dependence_group": g,
                "effect_size": float(row["effect_size"]),
                "standard_error": float(row["standard_error"]),
            })

    with tempfile.TemporaryDirectory() as tmpdir:
        in_tsv = Path(tmpdir) / "crosscheck_in.tsv"
        out_tsv = Path(tmpdir) / "crosscheck_r_out.tsv"
        pd.DataFrame(r_input_rows).to_csv(in_tsv, sep="\t", index=False)
        subprocess.run(
            ["Rscript", str(ROOT / "scripts" / "crosscheck_meta_analysis.R"), str(in_tsv), str(out_tsv)],
            check=True,
        )
        r_out = pd.read_csv(out_tsv, sep="\t")

    records = []
    max_abs = {m: 0.0 for m in ("pooled_effect", "pooled_standard_error", "tau2", "q", "i2_percent")}
    max_rel_p = 0.0
    max_mv_diff = 0.0

    for row in r_out.itertuples(index=False):
        cid = row.case_id
        method, rho, py = py_results[cid]
        diffs = {
            "pooled_effect": abs(py["pooled_effect"] - row.r_pooled_effect),
            "pooled_standard_error": abs(py["pooled_standard_error"] - row.r_pooled_se),
            "tau2": abs(py["tau2"] - row.r_tau2),
            "q": abs(py["q"] - row.r_q),
            "i2_percent": abs(py["i2_percent"] - row.r_i2_percent),
        }
        for k, v in diffs.items():
            if v > max_abs[k]:
                max_abs[k] = v
        rel_p = 0.0 if row.r_p_value == 0 and py["p_value"] == 0 else abs(py["p_value"] - row.r_p_value) / max(abs(row.r_p_value), 1e-300)
        if rel_p > max_rel_p:
            max_rel_p = rel_p

        mv_diff = float("nan")
        if not pd.isna(row.r_mv_pooled_effect):
            mv_diff = max(
                abs(py["pooled_effect"] - row.r_mv_pooled_effect),
                abs(py["pooled_standard_error"] - row.r_mv_pooled_se),
                abs(py["tau2"] - row.r_mv_tau2),
            )
            if mv_diff > max_mv_diff:
                max_mv_diff = mv_diff

        records.append({
            "case_id": cid,
            "method": method,
            "rho": "" if rho is None else rho,
            "py_pooled_effect": py["pooled_effect"],
            "r_pooled_effect": row.r_pooled_effect,
            "py_pooled_se": py["pooled_standard_error"],
            "r_pooled_se": row.r_pooled_se,
            "py_tau2": py["tau2"],
            "r_tau2": row.r_tau2,
            "py_i2": py["i2_percent"],
            "r_i2": row.r_i2_percent,
            "py_q": py["q"],
            "r_q": row.r_q,
            "py_p_value": py["p_value"],
            "r_p_value": row.r_p_value,
            "r_mv_pooled_effect": row.r_mv_pooled_effect,
            "r_mv_pooled_se": row.r_mv_pooled_se,
            "r_mv_tau2": row.r_mv_tau2,
            "max_abs_diff": max(diffs.values()),
            "rel_p_diff": rel_p,
            "mv_max_abs_diff": mv_diff,
        })

    out_table = pd.DataFrame(records)
    out_table.to_csv(output_dir / "metafor_crosscheck.tsv", sep="\t", index=False)

    # Assert strict numerical parity
    for metric, err in max_abs.items():
        tol = 1e-4 if metric == "i2_percent" else 1e-6
        if err > tol:
            raise AssertionError("Discrepancy {} exceeded tolerance {}: {}".format(metric, tol, err))
    if max_rel_p > 1e-5:
        raise AssertionError("Relative p-value discrepancy exceeded 1e-5: {}".format(max_rel_p))
    if max_mv_diff > 1e-5:
        raise AssertionError("rma.mv discrepancy exceeded 1e-5: {}".format(max_mv_diff))

    md = [
        "# R `metafor` Cross-Check Report",
        "",
        "Verified Python `aree.meta_analysis.random_effects` (DL, REML, and block-covariance shared-control GLS) against R `metafor` (`rma.uni` and `rma.mv`).",
        "",
        "- Cases tested: {} (12 synthetic edge/covariance cases + {} real-study gene groups from `data/harmonized/`).".format(
            len(out_table), len(out_table) - 12
        ),
        "- Max absolute difference in `pooled_effect`: `{:.3e}`".format(max_abs["pooled_effect"]),
        "- Max absolute difference in `pooled_standard_error`: `{:.3e}`".format(max_abs["pooled_standard_error"]),
        "- Max absolute difference in `tau2`: `{:.3e}`".format(max_abs["tau2"]),
        "- Max absolute difference in `q`: `{:.3e}`".format(max_abs["q"]),
        "- Max absolute difference in `i2_percent`: `{:.3e}`".format(max_abs["i2_percent"]),
        "- Max relative difference in `p_value`: `{:.3e}`".format(max_rel_p),
        "- Max absolute difference vs `metafor::rma.mv(yi, V, random = ~ 1 | group, method = 'REML')`: `{:.3e}`".format(max_mv_diff),
        "",
        "All checks passed within tolerance.",
        "",
    ]
    (output_dir / "metafor_crosscheck.md").write_text("\n".join(md))
    print("\n".join(md))
    return out_table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reports" / "synthesis")
    args = parser.parse_args()
    run_crosscheck(args.output_dir)
