from pathlib import Path

import pandas as pd

from aree.io import read_tsv
from aree.meta_analysis.random_effects import poolable
from aree.paths import root_path
from aree.prioritize.scoring import score_table
from aree.reporting.tables import dataframe_to_markdown


CARD_HEADER = "# Evidence Card: "
_WINDOWS_UNSAFE = str.maketrans({character: "_" for character in '<>:"/\\|?*'})
_WINDOWS_RESERVED = {"CON", "PRN", "AUX", "NUL"} | {"COM{}".format(i) for i in range(1, 10)} | {
    "LPT{}".format(i) for i in range(1, 10)
}


def card_filename(candidate_id):
    """Filename for a candidate's card that is valid on Windows, macOS and Linux."""
    name = str(candidate_id).translate(_WINDOWS_UNSAFE)
    name = "".join("_" if ord(character) < 32 else character for character in name).rstrip(" .") or "_"
    if name.split(".")[0].upper() in _WINDOWS_RESERVED:
        name = "_" + name
    return name + ".md"


def _check_filename_collisions(candidate_ids):
    # Checked before writing anything. Case-insensitive, because macOS and Windows filesystems
    # treat ABC.md and abc.md as the same file.
    claimed = {}
    for candidate_id in candidate_ids:
        key = card_filename(candidate_id).lower()
        if key in claimed:
            raise ValueError(
                "Candidates {!r} and {!r} map to the same card file {}".format(claimed[key], candidate_id, key)
            )
        claimed[key] = candidate_id


def _remove_stale_cards(output_dir, keep):
    # Only files this module generated (recognized by their header) are removed; anything else is left alone.
    for path in output_dir.glob("*.md"):
        if path not in keep and path.read_text().startswith(CARD_HEADER):
            path.unlink()


def build_evidence_cards(phenotype=None, evidence_path=None, scores_path=None, output_dir=None):
    """Write one card per feature.

    Without ``scores_path`` the cards are scored in memory from exactly the (filtered) evidence
    they display, so a phenotype filter never leaves a card showing a score from other contexts.
    Previously generated cards for candidates not in this build are removed, so the directory
    always matches the evidence it was built from.
    """
    evidence_path = evidence_path or root_path("data", "demo", "harmonized_evidence.tsv")
    output_dir = Path(output_dir) if output_dir else root_path("reports", "evidence_cards")
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence = read_tsv(evidence_path)
    if phenotype:
        evidence = evidence[evidence["phenotype"] == phenotype]
    scores = read_tsv(scores_path) if scores_path else score_table(evidence)
    _check_filename_collisions(evidence["feature_id_standardized"].unique())
    written = []
    for candidate_id, group in evidence.groupby("feature_id_standardized"):
        score_row = scores[scores["candidate_id"] == candidate_id]
        score_text = "not scored"
        category = "not scored"
        if not score_row.empty:
            score_text = str(score_row.iloc[0]["score"])
            category = score_row.iloc[0]["category"]
        # Explicit sort: value_counts orders ties differently across pandas versions.
        counts = group["molecular_direction"].value_counts()
        directions = ", ".join(
            "{}: {}".format(name, count) for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        )
        studies = ", ".join(sorted(group["study_id"].unique()))
        assays = ", ".join(sorted(group["feature_type"].unique()))
        contexts = "; ".join(
            sorted(set(group["phenotype"] + " / " + group["stressor"] + " / " + group["tissue"]))
        )
        limitations = "; ".join(sorted(set(group["quality_flags"].astype(str))))
        effect_summary = group[["study_id", "feature_type", "effect_size", "standard_error", "adjusted_p_value"]]
        body = [
            CARD_HEADER + str(candidate_id),
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
            "- Effects without standard errors (not pooled in meta-analysis): {}".format(
                ", ".join(sorted(group.loc[~poolable(group), "study_id"].unique())) or "none"
            ),
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
        path = output_dir / card_filename(candidate_id)
        path.write_text("\n".join(body))
        written.append(path)
    _remove_stale_cards(output_dir, set(written))
    return written
