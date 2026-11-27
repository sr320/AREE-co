"""Comparability-reviewed real cross-study evidence synthesis."""

import math
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import yaml

from aree.groups import column_groups, dependence_groups
from aree.io import read_evidence, read_tsv
from aree.meta_analysis.random_effects import EXPLORATORY_FLAG, RESULT_COLUMNS, _pool, inference_eligible
from aree.paths import root_path
from aree.reporting.tables import dataframe_to_markdown
from aree.validation.schemas import validate_study_file


SYNTHESIS_COLUMNS = [
    "stratum_id",
    "stratum_title",
    "phenotype",
    "stressor",
    "feature_id_standardized",
    "description",
    "feature_type",
    "effect_size_type",
    "method",
    "shared_control_rho",
    *RESULT_COLUMNS,
    "adjusted_p_value",
    "direction_concordant",
]


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR adjustment preserving missing values."""
    indexed = [(i, float(p)) for i, p in enumerate(p_values) if p is not None and not pd.isna(p)]
    adjusted = [float("nan")] * len(p_values)
    if not indexed:
        return adjusted
    indexed.sort(key=lambda item: item[1])
    n = len(indexed)
    running = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, pval = indexed[rank_idx]
        rank = rank_idx + 1
        running = min(running, pval * n / rank)
        adjusted[orig_idx] = min(1.0, running)
    return adjusted


def _load_studies(studies_dir=None):
    studies_dir = Path(studies_dir) if studies_dir else root_path("registry", "studies")
    studies = {}
    for path in sorted(studies_dir.glob("*.yaml")):
        study = validate_study_file(path)
        studies[study["study_id"]] = study
    return studies


def _load_annotations(analysis_dir=None):
    analysis_dir = Path(analysis_dir) if analysis_dir else root_path("reports", "analysis")
    descriptions = {}
    if not analysis_dir.exists():
        return descriptions
    for table_path in sorted(analysis_dir.glob("*/refseq_gene_annotations.tsv")):
        table = read_tsv(table_path, usecols=["feature_id_standardized", "description"])
        for feature, desc in zip(table["feature_id_standardized"], table["description"]):
            if feature not in descriptions and pd.notna(desc) and str(desc).strip():
                descriptions[feature] = str(desc)
    return descriptions


def load_synthesis_config(config_path=None, studies_dir=None, harmonized_dir=None):
    """Load and validate the comparability-reviewed synthesis configuration."""
    config_path = Path(config_path) if config_path else root_path("config", "synthesis_strata.yaml")
    config = yaml.safe_load(config_path.read_text())
    studies = _load_studies(studies_dir)

    seen_pooled = set()
    for stratum in config.get("strata", []):
        study_ids = stratum["study_ids"]
        if len(study_ids) < 2:
            raise ValueError("Stratum {} must include at least 2 contrasts".format(stratum["id"]))
        bioprojects = set()
        for sid in study_ids:
            if sid not in studies:
                raise ValueError("Unknown study {} in stratum {}".format(sid, stratum["id"]))
            study = studies[sid]
            if study["data_availability"]["status"].startswith("simulated"):
                raise ValueError("Simulated study {} cannot enter real synthesis".format(sid))
            if EXPLORATORY_FLAG in str(study.get("analysis_status", "")):
                raise ValueError("Exploratory study {} cannot enter real synthesis".format(sid))
            if sid in seen_pooled:
                raise ValueError("Study {} appears in multiple synthesis strata".format(sid))
            seen_pooled.add(sid)
            bioprojects.add((study.get("accessions") or {}).get("bioproject"))
        if len(bioprojects) < 2:
            raise ValueError("Stratum {} must span at least 2 independent BioProjects".format(stratum["id"]))

    unpooled = config.get("unpooled_contrasts", [])
    seen_unpooled = set()
    for item in unpooled:
        sid = item["study_id"]
        if sid not in studies:
            raise ValueError("Unknown unpooled study {}".format(sid))
        if sid in seen_pooled:
            raise ValueError("Study {} listed as both pooled and unpooled".format(sid))
        if sid in seen_unpooled:
            raise ValueError("Duplicate unpooled study {}".format(sid))
        seen_unpooled.add(sid)

    harmonized_dir = Path(harmonized_dir) if harmonized_dir else root_path("data", "harmonized")
    if harmonized_dir.exists():
        harmonized_ids = {p.stem for p in harmonized_dir.glob("*.tsv")}
        accounted = seen_pooled | seen_unpooled
        missing = sorted(harmonized_ids - accounted)
        if missing:
            raise ValueError(
                "Harmonized real studies missing from synthesis_strata.yaml: {}".format(", ".join(missing))
            )
    return config, studies


def synthesize_stratum(
    evidence,
    stratum,
    method="reml",
    shared_control_rho=0.5,
    annotations: Optional[Dict[str, str]] = None,
):
    """Run covariance-aware random-effects pooling within one comparability-reviewed stratum."""
    study_ids = list(stratum["study_ids"])
    sub = evidence[evidence["study_id"].isin(study_ids)].copy()
    present_studies = set(sub["study_id"].unique())
    if present_studies != set(study_ids):
        missing = sorted(set(study_ids) - present_studies)
        raise ValueError("Stratum {} missing evidence for studies: {}".format(stratum["id"], ", ".join(missing)))
    if not inference_eligible(sub).all():
        raise ValueError("Stratum {} contains non-inference-eligible exploratory rows".format(stratum["id"]))
    scales = set(sub["effect_size_type"].unique())
    if scales != {stratum["effect_size_type"]}:
        raise ValueError("Stratum {} has mismatched effect_size_type: {}".format(stratum["id"], sorted(scales)))

    sub["_inference_eligible"] = inference_eligible(sub)
    sub["_replication_group"] = dependence_groups(sub)
    n_indep_groups_stratum = sub["_replication_group"].nunique()
    if n_indep_groups_stratum < 2:
        raise ValueError("Stratum {} has fewer than 2 independent replication groups".format(stratum["id"]))

    annotations = annotations or {}
    pooled_rows = []
    total_features = 0
    single_group_features = 0
    groups = column_groups(
        sub,
        ["feature_id_standardized", "feature_type", "effect_size_type"],
        ["effect_size", "standard_error", "study_id", "_inference_eligible", "_replication_group"],
    )
    for (feature, feature_type, scale), columns in groups:
        total_features += 1
        res = _pool(
            columns["effect_size"],
            columns["standard_error"],
            columns["study_id"],
            columns["_inference_eligible"],
            columns["_replication_group"],
            method=method,
            shared_control_rho=shared_control_rho,
        )
        if res["n_independent_groups"] < 2:
            single_group_features += 1
            continue
        row = {
            "stratum_id": stratum["id"],
            "stratum_title": stratum["title"],
            "phenotype": stratum["phenotype"],
            "stressor": stratum["stressor"],
            "feature_id_standardized": feature,
            "description": annotations.get(feature, "Annotation unavailable"),
            "feature_type": feature_type,
            "effect_size_type": scale,
            "method": method,
            "shared_control_rho": shared_control_rho,
            **res,
        }
        pooled_rows.append(row)

    frame = pd.DataFrame(pooled_rows)
    if frame.empty:
        frame = pd.DataFrame(columns=SYNTHESIS_COLUMNS)
    else:
        frame["adjusted_p_value"] = benjamini_hochberg(frame["p_value"].tolist())
        frame["direction_concordant"] = frame["direction_consistency"] == 1.0
        frame["_abs_effect"] = frame["pooled_effect"].abs()
        frame = frame.sort_values(
            ["adjusted_p_value", "_abs_effect", "feature_id_standardized"],
            ascending=[True, False, True],
            na_position="last",
        ).drop(columns=["_abs_effect"])
        frame = frame[SYNTHESIS_COLUMNS].reset_index(drop=True)

    coverage = {
        "total_features_seen": total_features,
        "multi_group_pooled_features": len(frame),
        "single_group_only_features": single_group_features,
    }
    return frame, coverage


def _sensitivity_metrics(primary_df, alt_df, method, rho):
    merged = primary_df[["feature_id_standardized", "pooled_effect", "pooled_standard_error", "adjusted_p_value"]].merge(
        alt_df[["feature_id_standardized", "pooled_effect", "pooled_standard_error", "adjusted_p_value"]],
        on="feature_id_standardized",
        suffixes=("_primary", "_alt"),
    )
    prim_sig = set(primary_df.loc[primary_df["adjusted_p_value"] < 0.05, "feature_id_standardized"])
    alt_sig = set(alt_df.loc[alt_df["adjusted_p_value"] < 0.05, "feature_id_standardized"])
    union = prim_sig | alt_sig
    jaccard = len(prim_sig & alt_sig) / len(union) if union else 1.0
    alt_concordant = int(((alt_df["adjusted_p_value"] < 0.05) & alt_df["direction_concordant"]).sum())
    return {
        "stratum_id": str(primary_df["stratum_id"].iloc[0]),
        "method": method,
        "shared_control_rho": rho,
        "n_pooled": len(alt_df),
        "n_fdr_lt_0_05": len(alt_sig),
        "n_concordant_fdr_lt_0_05": alt_concordant,
        "median_i2_percent": round(float(alt_df["i2_percent"].median()), 2) if len(alt_df) else float("nan"),
        "median_se": round(float(alt_df["pooled_standard_error"].median()), 4) if len(alt_df) else float("nan"),
        "pearson_r_vs_primary": round(
            float(merged["pooled_effect_primary"].corr(merged["pooled_effect_alt"])), 5
        ) if len(merged) > 1 else float("nan"),
        "jaccard_fdr_0_05_vs_primary": round(jaccard, 4),
        "median_se_ratio_vs_primary": round(
            float((merged["pooled_standard_error_alt"] / merged["pooled_standard_error_primary"]).median()), 4
        ) if len(merged) else float("nan"),
    }


def _cross_context_heat_vs_selection(acute_heat_df, evidence):
    sel = evidence[evidence["study_id"] == "CGIG_THERMOTOL_RNASEQ_PRJNA694496"][
        ["feature_id_standardized", "effect_size", "standard_error", "adjusted_p_value"]
    ].rename(columns={
        "effect_size": "selection_log2fc",
        "standard_error": "selection_se",
        "adjusted_p_value": "selection_padj",
    })
    merged = acute_heat_df[[
        "feature_id_standardized", "description", "pooled_effect", "pooled_standard_error",
        "adjusted_p_value", "direction_concordant", "i2_percent",
    ]].merge(sel, on="feature_id_standardized", how="inner")
    valid = merged.dropna(subset=["pooled_effect", "selection_log2fc"])
    pearson = float(valid["pooled_effect"].corr(valid["selection_log2fc"], method="pearson"))
    spearman = float(valid["pooled_effect"].rank().corr(valid["selection_log2fc"].rank(), method="pearson"))
    both_sig = merged[(merged["adjusted_p_value"] < 0.05) & (merged["selection_padj"] < 0.05)]
    same_dir = both_sig[(both_sig["pooled_effect"] * both_sig["selection_log2fc"]) > 0]
    opp_dir = both_sig[(both_sig["pooled_effect"] * both_sig["selection_log2fc"]) < 0]
    return {
        "n_common": len(merged),
        "pearson_r": round(pearson, 4),
        "spearman_r": round(spearman, 4),
        "both_fdr_lt_0_05": len(both_sig),
        "both_sig_same_direction": len(same_dir),
        "both_sig_opposite_direction": len(opp_dir),
        "top_concordant": same_dir.sort_values("adjusted_p_value").head(15),
    }


def run_real_synthesis(
    evidence_path=None,
    config_path=None,
    data_dir=None,
    reports_dir=None,
    method=None,
    shared_control_rho=None,
):
    """Run comparability-reviewed synthesis across all curated strata and write reports + tables."""
    evidence_path = Path(evidence_path) if evidence_path else root_path("data", "harmonized")
    data_dir = Path(data_dir) if data_dir else root_path("data", "synthesis")
    reports_dir = Path(reports_dir) if reports_dir else root_path("reports", "synthesis")
    data_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    config, studies = load_synthesis_config(config_path=config_path, harmonized_dir=evidence_path if evidence_path.is_dir() else None)
    method = method or config.get("default_method", "reml")
    if shared_control_rho is None:
        shared_control_rho = float(config.get("default_shared_control_rho", 0.5))
    rhos = [float(r) for r in config.get("sensitivity_rhos", [0.0, 0.5, 0.8])]

    evidence = read_evidence(evidence_path)
    annotations = _load_annotations()

    stratum_summaries = []
    sensitivity_rows = []
    top_tables = []
    stratum_frames = {}

    for stratum in config["strata"]:
        primary_df, coverage = synthesize_stratum(
            evidence, stratum, method=method, shared_control_rho=shared_control_rho, annotations=annotations
        )
        stratum_frames[stratum["id"]] = primary_df
        out_tsv = data_dir / "{}_meta.tsv".format(stratum["id"])
        primary_df.to_csv(out_tsv, sep="\t", index=False)

        sub_ev = evidence[evidence["study_id"].isin(stratum["study_ids"])]
        dep_groups = set(dependence_groups(sub_ev))
        bioprojects = {(studies[sid].get("accessions") or {}).get("bioproject") for sid in stratum["study_ids"]}
        sig = primary_df[primary_df["adjusted_p_value"] < 0.05]
        concordant_sig = sig[sig["direction_concordant"]]

        stratum_summaries.append({
            "stratum_id": stratum["id"],
            "title": stratum["title"],
            "phenotype": stratum["phenotype"],
            "stressor": stratum["stressor"],
            "n_contrasts": len(stratum["study_ids"]),
            "n_bioprojects": len(bioprojects),
            "n_independent_groups": len(dep_groups),
            "total_features_seen": coverage["total_features_seen"],
            "multi_group_pooled_features": coverage["multi_group_pooled_features"],
            "single_group_only_features": coverage["single_group_only_features"],
            "pooled_fdr_lt_0_05": len(sig),
            "concordant_fdr_lt_0_05": len(concordant_sig),
            "concordant_positive": int((concordant_sig["pooled_effect"] > 0).sum()),
            "concordant_negative": int((concordant_sig["pooled_effect"] < 0).sum()),
            "median_i2_percent": round(float(primary_df["i2_percent"].median()), 2) if len(primary_df) else float("nan"),
            "study_ids": ";".join(stratum["study_ids"]),
        })

        # Sensitivity grid: (method, shared_control_rho), plus DL at primary rho
        configs_to_check = [(method, shared_control_rho)]
        alt_method = "dl" if method == "reml" else "reml"
        configs_to_check.append((alt_method, shared_control_rho))
        for alt_rho in rhos:
            if (method, alt_rho) not in configs_to_check:
                configs_to_check.append((method, alt_rho))
        for alt_m, alt_r in configs_to_check:
            if alt_m == method and math.isclose(alt_r, shared_control_rho):
                sensitivity_rows.append(_sensitivity_metrics(primary_df, primary_df, alt_m, alt_r))
            else:
                alt_df, _ = synthesize_stratum(
                    evidence, stratum, method=alt_m, shared_control_rho=alt_r, annotations=annotations
                )
                sensitivity_rows.append(_sensitivity_metrics(primary_df, alt_df, alt_m, alt_r))

        # Bounded top 25 pooled associations per stratum (prioritizing concordant FDR < 0.05 first)
        top = primary_df.assign(_concordant=primary_df["direction_concordant"].astype(int)).sort_values(
            ["_concordant", "adjusted_p_value", "feature_id_standardized"],
            ascending=[False, True, True],
        ).drop(columns=["_concordant"]).head(25)
        top_tables.append(top)

    summary_df = pd.DataFrame(stratum_summaries)
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    top_df = pd.concat(top_tables, ignore_index=True) if top_tables else pd.DataFrame(columns=SYNTHESIS_COLUMNS)

    summary_df.to_csv(reports_dir / "stratum_summary.tsv", sep="\t", index=False)
    sensitivity_df.to_csv(reports_dir / "sensitivity_summary.tsv", sep="\t", index=False)
    top_df.to_csv(reports_dir / "top_pooled_associations.tsv", sep="\t", index=False)

    # Cross-context non-pooled comparison: acute gill heat shock vs PRJNA694496 larval selection
    heat_cmp = None
    if "acute_heat_response_gill" in stratum_frames and "CGIG_THERMOTOL_RNASEQ_PRJNA694496" in set(evidence["study_id"]):
        heat_cmp = _cross_context_heat_vs_selection(stratum_frames["acute_heat_response_gill"], evidence)
        heat_cmp["top_concordant"].to_csv(reports_dir / "acute_heat_vs_prjna694496_concordant.tsv", sep="\t", index=False)

    md_parts = [
        "# Real-Study Cross-Study Synthesis Report",
        "",
        "**Association evidence only; not validated biomarkers.** Cross-study pooling is restricted to comparability-reviewed biological strata with at least two independent BioProjects and two independent replication groups, using `{}` random-effects meta-analysis and shared-control block covariance (`rho = {}`).".format(
            method.upper(), shared_control_rho
        ),
        "",
        "## Stratum Summary",
        "",
        dataframe_to_markdown(summary_df[[
            "stratum_id", "n_contrasts", "n_bioprojects", "n_independent_groups",
            "multi_group_pooled_features", "pooled_fdr_lt_0_05", "concordant_fdr_lt_0_05", "median_i2_percent",
        ]]),
        "",
        "## Estimator and Covariance Sensitivity",
        "",
        dataframe_to_markdown(sensitivity_df[[
            "stratum_id", "method", "shared_control_rho", "n_fdr_lt_0_05",
            "n_concordant_fdr_lt_0_05", "median_i2_percent", "pearson_r_vs_primary", "jaccard_fdr_0_05_vs_primary",
        ]]),
        "",
    ]

    for stratum in config["strata"]:
        sid = stratum["id"]
        sub_top = top_df[top_df["stratum_id"] == sid].head(10)
        md_parts.extend([
            "## {} (`{}`)".format(stratum["title"], sid),
            "",
            stratum["comparability_rationale"],
            "",
            "**Limitations:**",
            *("- " + lim for lim in stratum["limitations"]),
            "",
            "### Top 10 pooled associations (concordant signals prioritized)",
            "",
            dataframe_to_markdown(sub_top[[
                "feature_id_standardized", "description", "n_studies", "n_independent_groups",
                "pooled_effect", "pooled_standard_error", "adjusted_p_value", "i2_percent", "direction_concordant",
            ]]),
            "",
        ])

    if heat_cmp is not None:
        md_parts.extend([
            "## Cross-Context Comparison: Acute Heat Gill Meta-Analysis vs PRJNA694496 Constitutive Selection (Not Pooled)",
            "",
            "Consistent with AREC's biological boundaries, constitutive larval expression in thermotolerance-selected oysters (`CGIG_THERMOTOL_RNASEQ_PRJNA694496`) is **not** pooled with acute adult/juvenile gill heat-shock responses (`acute_heat_response_gill`). Across {:,} common genes: Pearson $r = {:.4f}$, Spearman $\\rho = {:.4f}$; {:,} genes reach FDR < 0.05 in both ({:,} same direction, {:,} opposite direction).".format(
                heat_cmp["n_common"], heat_cmp["pearson_r"], heat_cmp["spearman_r"],
                heat_cmp["both_fdr_lt_0_05"], heat_cmp["both_sig_same_direction"], heat_cmp["both_sig_opposite_direction"],
            ),
            "",
        ])

    md_parts.extend([
        "## Unpooled Harmonized Contrasts",
        "",
        "Ten harmonized contrasts are intentionally withheld from cross-study pooling:",
        "",
        *("- `{}`: {}".format(item["study_id"], item["reason"]) for item in config["unpooled_contrasts"]),
        "",
    ])

    (reports_dir / "synthesis_report.md").write_text("\n".join(md_parts))
    return reports_dir / "synthesis_report.md"
