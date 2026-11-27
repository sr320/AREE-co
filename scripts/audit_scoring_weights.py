#!/usr/bin/env python3
"""Audit candidate scoring weights and quality-flag classification across real-study strata."""

import argparse
from pathlib import Path

import pandas as pd

from aree.io import read_evidence
from aree.meta_analysis.synthesis import load_synthesis_config
from aree.prioritize import scoring
from aree.reporting.tables import dataframe_to_markdown

ROOT = Path(__file__).resolve().parents[1]


def _score_with_weights(evidence, weights):
    original = dict(scoring.SCORE_WEIGHTS)
    original_total = scoring._WEIGHT_TOTAL
    try:
        scoring.SCORE_WEIGHTS.clear()
        scoring.SCORE_WEIGHTS.update(weights)
        scoring._WEIGHT_TOTAL = sum(weights.values())
        return scoring.score_table(evidence, method="reml", shared_control_rho=0.5)
    finally:
        scoring.SCORE_WEIGHTS.clear()
        scoring.SCORE_WEIGHTS.update(original)
        scoring._WEIGHT_TOTAL = original_total


def _compare_rankings(base_df, alt_df, stratum_id, scheme):
    merged = base_df[["candidate_id", "score"]].merge(
        alt_df[["candidate_id", "score"]], on="candidate_id", suffixes=("_base", "_alt")
    )
    spearman = float(merged["score_base"].rank().corr(merged["score_alt"].rank(), method="pearson"))
    top25_base = set(base_df.head(25)["candidate_id"])
    top25_alt = set(alt_df.head(25)["candidate_id"])
    top100_base = set(base_df.head(100)["candidate_id"])
    top100_alt = set(alt_df.head(100)["candidate_id"])
    j25 = len(top25_base & top25_alt) / len(top25_base | top25_alt) if (top25_base | top25_alt) else 1.0
    j100 = len(top100_base & top100_alt) / len(top100_base | top100_alt) if (top100_base | top100_alt) else 1.0
    return {
        "stratum_id": stratum_id,
        "scheme": scheme,
        "n_candidates": len(base_df),
        "spearman_r_vs_default": round(spearman, 4),
        "top25_jaccard": round(j25, 4),
        "top100_jaccard": round(j100, 4),
    }


def run_audit(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config, _ = load_synthesis_config()
    evidence = read_evidence(ROOT / "data" / "harmonized")

    # Audit quality flags across all real evidence
    flag_rows = []
    all_flags = sorted({
        tok.strip()
        for flags in evidence["quality_flags"].dropna().unique()
        for tok in str(flags).split(";")
        if tok.strip()
    })
    for tok in all_flags:
        flag_rows.append({
            "quality_flag": tok,
            "classified_as": "provenance_non_penalty" if scoring._is_non_penalty_flag(tok) else "caveat_penalty",
        })
    flag_df = pd.DataFrame(flag_rows)

    default_weights = dict(scoring.SCORE_WEIGHTS)
    schemes = {
        "equal_weights": {k: 1.0 for k in default_weights},
        "replication_heavy_2x": {
            **default_weights,
            "n_studies": default_weights["n_studies"] * 2.0,
            "direction_consistency": default_weights["direction_consistency"] * 2.0,
        },
        "outcome_quality_heavy_2x": {
            **default_weights,
            "phenotype_relevance": default_weights["phenotype_relevance"] * 2.0,
            "data_quality": default_weights["data_quality"] * 2.0,
        },
    }
    for key in default_weights:
        schemes["drop_" + key] = {k: (0.0 if k == key else v) for k, v in default_weights.items()}

    audit_rows = []
    for stratum in config["strata"]:
        sub = evidence[evidence["study_id"].isin(stratum["study_ids"])].copy()
        base_df = scoring.score_table(sub, method="reml", shared_control_rho=0.5)
        for name, weights in schemes.items():
            alt_df = _score_with_weights(sub, weights)
            audit_rows.append(_compare_rankings(base_df, alt_df, stratum["id"], name))

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(output_dir / "scoring_weight_audit.tsv", sep="\t", index=False)
    flag_df.to_csv(output_dir / "quality_flag_classification.tsv", sep="\t", index=False)

    summary_by_scheme = audit_df.groupby("scheme", sort=False).agg(
        mean_spearman=("spearman_r_vs_default", "mean"),
        min_spearman=("spearman_r_vs_default", "min"),
        mean_top25_jaccard=("top25_jaccard", "mean"),
        mean_top100_jaccard=("top100_jaccard", "mean"),
    ).reset_index()

    md = [
        "# Candidate Scoring Weight & Quality-Flag Audit",
        "",
        "Evaluated `score_table(..., method='reml', shared_control_rho=0.5)` across all {} comparability-reviewed strata under equal weights, replication-heavy, outcome/quality-heavy, and leave-one-component-out perturbations.".format(
            len(config["strata"])
        ),
        "",
        "## Weight Stability Summary",
        "",
        dataframe_to_markdown(summary_by_scheme.round(4)),
        "",
        "## Quality-Flag Classification Across All Harmonized Real Studies",
        "",
        dataframe_to_markdown(flag_df),
        "",
    ]
    (output_dir / "scoring_weight_audit.md").write_text("\n".join(md))
    print("Wrote scoring audit to", output_dir)
    return audit_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reports" / "synthesis")
    args = parser.parse_args()
    run_audit(args.output_dir)
