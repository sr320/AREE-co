from pathlib import Path

import pandas as pd
import streamlit as st

from aree.io import read_tsv
from aree.reporting.evidence_cards import card_filename
from aree.reporting.site_findings import load_findings, number
from aree.reporting.tables import dataframe_to_markdown
from aree.validation.schemas import validate_study_file


ROOT = Path(__file__).resolve().parents[1]
# Streamlit 1.49 replaced use_container_width (since deprecated) with width="stretch"; Python 3.8
# tops out at older Streamlit releases, so pick whichever keyword the installed version supports.
_STREAMLIT_VERSION = tuple(int(part) for part in st.__version__.split(".")[:2])
FULL_WIDTH = {"width": "stretch"} if _STREAMLIT_VERSION >= (1, 49) else {"use_container_width": True}
_cache_data = getattr(st, "cache_data", lambda fn: fn)
MAX_DISPLAY_ROWS = 500
REAL_EVIDENCE_DISPLAY_COLUMNS = [
    "feature_id_standardized",
    "description",
    "study_id",
    "sample_comparison",
    "tissue",
    "life_stage",
    "stressor",
    "effect_size_type",
    "effect_size",
    "standard_error",
    "adjusted_p_value",
    "inference_status",
    "quality_flags",
]
REAL_SEARCH_COLUMNS = [
    "feature_id_standardized",
    "description",
    "study_id",
    "sample_comparison",
    "quality_flags",
]
HIGHLIGHT_COLUMNS = [
    "study_id",
    "feature_id_standardized",
    "description",
    "sample_comparison",
    "effect_size_type",
    "effect_size",
    "standard_error",
    "adjusted_p_value",
    "inference_status",
    "selection_basis",
]


@_cache_data
def _load_registry(registry_path):
    return pd.read_csv(registry_path)


@_cache_data
def _load_tsv(path_str):
    return read_tsv(path_str)


@_cache_data
def _load_real_bundle(root_str):
    root = Path(root_str)
    studies = {}
    for path in sorted((root / "registry" / "studies").glob("*.yaml")):
        study = validate_study_file(path)
        studies[study["study_id"]] = study
    evidence, summaries, highlights = load_findings(root, studies)
    return studies, evidence, summaries, highlights


def _format_real_card(feature, rows, studies):
    parts = [
        "# Real Evidence Card: {}".format(feature),
        "",
        "**Real public-study evidence. Association evidence only; not a validated biomarker.**",
        "",
        (
            "Differential expression alone does not establish a resilience marker; higher expression "
            "does not necessarily mean greater resilience. All available harmonized study rows for "
            "this feature are shown below, including nonsignificant and exploratory results. "
            "No cross-study score or pooled effect is calculated."
        ),
        "",
        "## Study-Specific Effects",
        "",
    ]
    effect_rows = []
    for row in rows.itertuples(index=False):
        study = studies[row.study_id]
        bioproject = (study.get("accessions") or {}).get("bioproject", "")
        exploratory = str(row.inference_status).startswith("Exploratory")
        adj_p = number(row.adjusted_p_value) + (" (nominal only)" if exploratory else "")
        effect_rows.append(
            {
                "Study": row.study_id,
                "BioProject": bioproject,
                "RefSeq description": row.description,
                "Contrast": row.sample_comparison,
                "Effect scale": row.effect_size_type.replace("_", " "),
                "Effect": number(row.effect_size),
                "Standard error": number(row.standard_error),
                "Adjusted p-value": adj_p,
                "Interpretation": row.inference_status,
            }
        )
    parts.append(dataframe_to_markdown(pd.DataFrame(effect_rows)))
    for row in rows.itertuples(index=False):
        study = studies[row.study_id]
        parts.extend(
            [
                "",
                "## {}".format(row.study_id),
                "",
                "- **Phenotype interpretation:** {}".format(study["phenotype_direction"]),
                "- **Biological replication:** {}".format(study["biological_replication"]),
                "- **Evidence classification:** {}. Tissue: {}; life stage: {}.".format(
                    row.resilience_classification, row.tissue, row.life_stage
                ),
                "- **Quality flags:** `{}`".format(row.quality_flags),
            ]
        )
        if str(row.inference_status).startswith("Exploratory"):
            parts.append(
                "- **Exploratory uncertainty:** Tank clustering is unmodeled. Standard errors and "
                "adjusted p-values above are nominal only; excluded from pooling and significance rewards."
            )
        parts.append("- **Study limitations:**")
        parts.extend("  - {}".format(text) for text in study["limitations"])
        parts.append(
            "- **Provenance:** Reference `{}` / `{}`; mapping `{}`; workflow `{}`; input SHA-256 `{}`.".format(
                row.genome_assembly, row.annotation_version, row.mapping_confidence,
                row.workflow_version, row.input_checksum,
            )
        )
    return "\n".join(parts)


st.set_page_config(page_title="AREE", layout="wide")
st.title("Aquaculture Resilience Evidence Engine")
mode = st.radio(
    "Dataset mode",
    ["Real public studies", "Synthetic demo"],
    horizontal=True,
)
real_mode = mode == "Real public studies"
if real_mode:
    st.info(
        "Real public-study evidence below presents study-level molecular associations only — "
        "not validated biomarkers. No cross-study pooled effect or global candidate score is computed; "
        "exploratory unmodeled-tank p-values are labeled nominal only."
    )
else:
    st.warning(
        "Evidence, candidate scores and cards below use simulated demo data. "
        "They are not findings from the registered real public studies."
    )

registry_path = ROOT / "registry" / "study_registry.csv"
demo_evidence_path = ROOT / "data" / "demo" / "harmonized_evidence.tsv"
demo_scores_path = ROOT / "data" / "demo" / "candidate_scores.tsv"

tabs = st.tabs(["Studies", "Evidence", "Findings & Cards", "Pipeline Status"])

with tabs[0]:
    st.subheader("Registered Studies")
    if registry_path.exists():
        registry = _load_registry(str(registry_path))
        is_sim = registry["data_status"].astype(str).str.startswith("simulated")
        filtered_registry = registry[~is_sim] if real_mode else registry[is_sim]
        st.dataframe(filtered_registry, **FULL_WIDTH)
        filename = "real_study_registry.csv" if real_mode else "demo_study_registry.csv"
        st.download_button("Download registry CSV", filtered_registry.to_csv(index=False), filename)
    else:
        st.info("No registry CSV found.")

with tabs[1]:
    st.subheader("Harmonized Molecular Evidence")
    if real_mode:
        studies, evidence, summaries, highlights = _load_real_bundle(str(ROOT))
        cols = st.columns(4)
        stressor = cols[0].multiselect("Stressor", sorted(evidence["stressor"].dropna().unique()))
        tissue = cols[1].multiselect("Tissue", sorted(evidence["tissue"].dropna().unique()))
        study_filter = cols[2].multiselect("Study contrast", sorted(evidence["study_id"].dropna().unique()))
        inference = cols[3].selectbox("Inference", ["All studies", "Inference-eligible only", "Exploratory only"])
        query = st.text_input("Search gene ID, RefSeq description, contrast, or quality flags")

        filtered = evidence
        for column, values in [("stressor", stressor), ("tissue", tissue), ("study_id", study_filter)]:
            if values:
                filtered = filtered[filtered[column].isin(values)]
        if inference == "Inference-eligible only":
            filtered = filtered[~filtered["inference_status"].astype(str).str.startswith("Exploratory")]
        elif inference == "Exploratory only":
            filtered = filtered[filtered["inference_status"].astype(str).str.startswith("Exploratory")]
        if query:
            mask = pd.Series(False, index=filtered.index)
            for column in REAL_SEARCH_COLUMNS:
                mask = mask | filtered[column].astype(str).str.contains(query, case=False, na=False, regex=False)
            filtered = filtered[mask]

        if inference == "Exploratory only":
            filtered = filtered.assign(_abs_effect=filtered["effect_size"].abs()).sort_values(
                ["_abs_effect", "feature_id_standardized"], ascending=[False, True], na_position="last"
            ).drop(columns=["_abs_effect"])
        else:
            filtered = filtered.assign(
                _exploratory=filtered["inference_status"].astype(str).str.startswith("Exploratory"),
                _abs_effect=filtered["effect_size"].abs(),
            ).sort_values(
                ["_exploratory", "adjusted_p_value", "_abs_effect", "feature_id_standardized"],
                ascending=[True, True, False, True],
                na_position="last",
            ).drop(columns=["_exploratory", "_abs_effect"])

        display = filtered[REAL_EVIDENCE_DISPLAY_COLUMNS]
        total_matches = len(display)
        if total_matches > MAX_DISPLAY_ROWS:
            st.caption(
                "{:,} matching records (showing first {:,}; download TSV for the full filtered set).".format(
                    total_matches, MAX_DISPLAY_ROWS
                )
            )
        else:
            st.caption("{:,} matching records.".format(total_matches))
        st.dataframe(display.head(MAX_DISPLAY_ROWS), **FULL_WIDTH)
        st.download_button(
            "Download filtered evidence TSV",
            display.to_csv(sep="\t", index=False),
            "aree_real_filtered_evidence.tsv",
        )
    elif demo_evidence_path.exists():
        evidence = _load_tsv(str(demo_evidence_path))
        cols = st.columns(5)
        stressor = cols[0].multiselect("Stressor", sorted(evidence["stressor"].dropna().unique()))
        phenotype = cols[1].multiselect("Phenotype", sorted(evidence["phenotype"].dropna().unique()))
        assay = cols[2].multiselect("Feature type", sorted(evidence["feature_type"].dropna().unique()))
        tissue = cols[3].multiselect("Tissue", sorted(evidence["tissue"].dropna().unique()))
        species = cols[4].multiselect("Species", sorted(evidence["species"].dropna().unique()))
        query = st.text_input("Search gene, protein, orthogroup, pathway, or original ID")
        filtered = evidence.copy()
        for column, values in [
            ("stressor", stressor), ("phenotype", phenotype), ("feature_type", assay),
            ("tissue", tissue), ("species", species),
        ]:
            if values:
                filtered = filtered[filtered[column].isin(values)]
        if query:
            mask = filtered.astype(str).apply(
                lambda col: col.str.contains(query, case=False, na=False, regex=False)
            ).any(axis=1)
            filtered = filtered[mask]
        st.dataframe(filtered, **FULL_WIDTH)
        st.download_button(
            "Download filtered evidence TSV",
            filtered.to_csv(sep="\t", index=False),
            "aree_filtered_evidence.tsv",
        )
    else:
        st.info("Run `aree harmonize-demo` to generate evidence.")

with tabs[2]:
    if real_mode:
        st.subheader("Real-Study Findings & Evidence Cards")
        st.caption(
            "Up to ten displayed findings per harmonized study for browsing; no global biomarker score "
            "or cross-study pooled effect is calculated."
        )
        st.dataframe(summaries, **FULL_WIDTH)
        st.dataframe(highlights[HIGHLIGHT_COLUMNS], **FULL_WIDTH)
        featured_genes = sorted(highlights["feature_id_standardized"].astype(str).unique())
        if featured_genes:
            selected = st.selectbox("Real evidence card", featured_genes)
            card_rows = evidence[evidence["feature_id_standardized"] == selected]
            st.markdown(_format_real_card(selected, card_rows, studies))
    else:
        st.subheader("Simulated Demo Candidate Scores & Evidence Cards")
        if demo_scores_path.exists():
            scores = _load_tsv(str(demo_scores_path))
            st.dataframe(scores, **FULL_WIDTH)
            selected = st.selectbox("Evidence card", sorted(scores["candidate_id"].astype(str)))
            card = ROOT / "reports" / "evidence_cards" / card_filename(selected)
            if card.exists():
                st.markdown(card.read_text())
        else:
            st.info("Run `aree build-evidence-cards` to generate candidates.")

with tabs[3]:
    st.subheader("Pipeline Status")
    status = pd.DataFrame(
        [
            {"layer": "Study registry & controlled vocabularies", "status": "complete and runnable"},
            {"layer": "Processed harmonization", "status": "complete and runnable for MVP schema"},
            {
                "layer": "Raw RNA-seq Nextflow pipeline",
                "status": "complete (paired-end, single-end, covariates, collapse-by, 3' tag-seq)",
            },
            {"layer": "Raw methylation/proteomics/metabolomics pipelines", "status": "planned; processed results only"},
            {"layer": "Meta-analysis and scoring (demo / configurable)", "status": "complete and runnable"},
            {"layer": "Static Quarto site & DuckDB-WASM Query page", "status": "complete and deployed via Pages"},
        ]
    )
    st.dataframe(status, **FULL_WIDTH)
