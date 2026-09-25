from pathlib import Path

import pandas as pd

from aree.io import read_tsv
from aree.paths import root_path
from aree.reporting.tables import dataframe_to_markdown


def build_demo_report(output_path=None, registry_path=None, evidence_path=None, scores_path=None):
    output_path = Path(output_path) if output_path else root_path("reports", "demo_report.md")
    registry = pd.read_csv(registry_path or root_path("registry", "study_registry.csv"))
    evidence = read_tsv(evidence_path or root_path("data", "demo", "harmonized_evidence.tsv"))
    scores = read_tsv(scores_path or root_path("data", "demo", "candidate_scores.tsv"))
    body = [
        "# AREE Demo Report",
        "",
        "This report uses simulated demo evidence to demonstrate the MVP evidence-generation workflow.",
        "",
        "## Registered Studies",
        "",
        dataframe_to_markdown(registry),
        "",
        "## Evidence Counts",
        "",
        dataframe_to_markdown(evidence.groupby(["phenotype", "stressor", "feature_type"]).size().reset_index(name="records")),
        "",
        "## Candidate Scores",
        "",
        dataframe_to_markdown(scores),
        "",
        "## Interpretation Guardrail",
        "",
        "Candidates are prioritized associations. A statistically significant single-study result is not treated as validation.",
        "",
    ]
    output_path.write_text("\n".join(body))
    return output_path
