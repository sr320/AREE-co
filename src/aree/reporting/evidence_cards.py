from collections import Counter
from pathlib import Path

from aree.groups import column_groups, is_missing, present
from aree.io import read_tsv
from aree.meta_analysis.random_effects import poolable
from aree.paths import root_path
from aree.prioritize.scoring import score_table
from aree.reporting.tables import rows_to_markdown


CARD_HEADER = "# Evidence Card: "
EFFECT_SUMMARY_COLUMNS = ["study_id", "feature_type", "effect_size", "standard_error", "adjusted_p_value"]
CARD_COLUMNS = EFFECT_SUMMARY_COLUMNS + [
    "molecular_direction",
    "phenotype",
    "stressor",
    "tissue",
    "species",
    "ortholog_reference",
    "mapping_confidence",
    "quality_flags",
    "_poolable",
]
_WINDOWS_UNSAFE = str.maketrans({character: "_" for character in '<>:"/\\|?*'})
_WINDOWS_RESERVED = {"CON", "PRN", "AUX", "NUL"} | {"COM{}".format(i) for i in range(1, 10)} | {
    "LPT{}".format(i) for i in range(1, 10)
}


def _joined(values):
    return ", ".join(sorted({str(value) for value in present(values)}))


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
    scores = scores.drop_duplicates("candidate_id")
    score_lookup = dict(zip(scores["candidate_id"], zip(scores["score"], scores["category"])))
    _check_filename_collisions(evidence["feature_id_standardized"].unique())
    usable = poolable(evidence).to_numpy()
    written = []
    groups = column_groups(evidence.assign(_poolable=usable), "feature_id_standardized", CARD_COLUMNS)
    for candidate_id, group in groups:
        score_text, category = "not scored", "not scored"
        if candidate_id in score_lookup:
            score, category = score_lookup[candidate_id]
            score_text = str(score)
        # Explicit sort: counts must not depend on how a library orders ties.
        counts = Counter(present(group["molecular_direction"]))
        directions = ", ".join(
            "{}: {}".format(name, count) for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        )
        contexts = "; ".join(
            sorted(
                {
                    "{} / {} / {}".format(phenotype, stressor, tissue)
                    for phenotype, stressor, tissue in zip(group["phenotype"], group["stressor"], group["tissue"])
                    if not any(is_missing(part) for part in (phenotype, stressor, tissue))
                }
            )
        )
        not_pooled = sorted({study for study, ok in zip(group["study_id"], group["_poolable"]) if not ok})
        effect_rows = list(zip(*(group[column] for column in EFFECT_SUMMARY_COLUMNS)))
        body = [
            CARD_HEADER + str(candidate_id),
            "",
            "**Status:** Association evidence only; not a validated biomarker.",
            "",
            "- Candidate score: {}".format(score_text),
            "- Ranking category: {}".format(category),
            "- Species context: {}".format(_joined(group["species"])),
            "- Ortholog/reference context: {}".format(
                ", ".join(sorted({str(value) for value in present(group["ortholog_reference"])})) or "not resolved"
            ),
            "- Supporting studies: {}".format(_joined(group["study_id"])),
            "- Assay types represented: {}".format(_joined(group["feature_type"])),
            "- Phenotype/stressor/tissue contexts: {}".format(contexts),
            "- Direction of association: {}".format(directions),
            "- Identifier mapping confidence: {}".format(_joined(group["mapping_confidence"])),
            "- Limitations: {}".format("; ".join(sorted({str(value) for value in group["quality_flags"]}))),
            "- Effects without standard errors (not pooled in meta-analysis): {}".format(
                ", ".join(not_pooled) or "none"
            ),
            "",
            "## Effect Summary",
            "",
            rows_to_markdown(EFFECT_SUMMARY_COLUMNS, effect_rows),
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
