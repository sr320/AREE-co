from pathlib import Path

import pandas as pd

from aree.paths import root_path
from aree.reporting.tables import dataframe_to_markdown


def build_evidence_cards(phenotype=None, evidence_path=None, scores_path=None, output_dir=None):
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    scores_path = scores_path or root_path("data", "demo", "candidate_scores.tsv")
    output_dir = Path(output_dir) if output_dir else root_path("reports", "evidence_cards")
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_csv(evidence_path, sep="\t")
    scores = pd.read_csv(scores_path, sep="\t")
    if phenotype:
        evidence = evidence[evidence["phenotype"] == phenotype]
    written = []
    for candidate_id, group in evidence.groupby("feature_id_standardized"):
        score_row = scores[scores["candidate_id"] == candidate_id]
        score_text = "not scored"
        category = "not scored"
        if not score_row.empty:
            score_text = str(score_row.iloc[0]["score"])
            category = score_row.iloc[0]["category"]
        directions = group["molecular_direction"].value_counts().to_dict()
        studies = ", ".join(sorted(group["study_id"].unique()))
        assays = ", ".join(sorted(group["feature_type"].unique()))
        contexts = "; ".join(
            sorted(set(group["phenotype"] + " / " + group["stressor"] + " / " + group["tissue"]))
        )
        limitations = "; ".join(sorted(set(group["quality_flags"].astype(str))))
        effect_summary = group[["study_id", "feature_type", "effect_size", "standard_error", "adjusted_p_value"]]
        body = [
            "# Evidence Card: {}".format(candidate_id),
            "",
            "**Status:** Association evidence only; not a validated biomarker.",
            "",
            "- Candidate score: {}".format(score_text),
            "- Ranking category: {}".format(category),
            "- Species context: {}".format(", ".join(sorted(group["species"].unique()))),
            "- Ortholog/reference context: {}".format(", ".join(sorted(set(group["ortholog_reference"].dropna().astype(str)))) or "not resolved"),
            "- Supporting studies: {}".format(studies),
            "- Assay types represented: {}".format(assays),
            "- Phenotype/stressor/tissue contexts: {}".format(contexts),
            "- Direction of association: {}".format(directions),
            "- Identifier mapping confidence: {}".format(", ".join(sorted(group["mapping_confidence"].unique()))),
            "- Limitations: {}".format(limitations),
            "",
            "## Effect Summary",
            "",
            dataframe_to_markdown(effect_summary),
            "",
            "## Recommended Next Validation Step",
            "",
            "Prioritize independent biological replication with matched phenotype definitions and targeted validation in the relevant tissue and life stage.",
            "",
        ]
        path = output_dir / "{}.md".format(candidate_id.replace("/", "_"))
        path.write_text("\n".join(body))
        written.append(path)
    return written
