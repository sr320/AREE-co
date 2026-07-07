from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parents[1]


st.set_page_config(page_title="AREE", layout="wide")
st.title("Aquaculture Resilience Evidence Engine")

registry_path = ROOT / "registry" / "study_registry.csv"
evidence_path = ROOT / "data" / "demo" / "harmonized_evidence.tsv"
scores_path = ROOT / "data" / "demo" / "candidate_scores.tsv"

tabs = st.tabs(["Studies", "Evidence", "Candidates", "Pipeline Status"])

with tabs[0]:
    st.subheader("Registered Studies")
    if registry_path.exists():
        registry = pd.read_csv(registry_path)
        st.dataframe(registry, use_container_width=True)
        st.download_button("Download registry CSV", registry.to_csv(index=False), "study_registry.csv")
    else:
        st.info("No registry CSV found.")

with tabs[1]:
    st.subheader("Harmonized Molecular Evidence")
    if evidence_path.exists():
        evidence = pd.read_csv(evidence_path, sep="\t")
        cols = st.columns(5)
        stressor = cols[0].multiselect("Stressor", sorted(evidence["stressor"].dropna().unique()))
        phenotype = cols[1].multiselect("Phenotype", sorted(evidence["phenotype"].dropna().unique()))
        assay = cols[2].multiselect("Feature type", sorted(evidence["feature_type"].dropna().unique()))
        tissue = cols[3].multiselect("Tissue", sorted(evidence["tissue"].dropna().unique()))
        species = cols[4].multiselect("Species", sorted(evidence["species"].dropna().unique()))
        query = st.text_input("Search gene, protein, orthogroup, pathway, or original ID")
        filtered = evidence.copy()
        for column, values in [("stressor", stressor), ("phenotype", phenotype), ("feature_type", assay), ("tissue", tissue), ("species", species)]:
            if values:
                filtered = filtered[filtered[column].isin(values)]
        if query:
            mask = filtered.astype(str).apply(lambda col: col.str.contains(query, case=False, na=False)).any(axis=1)
            filtered = filtered[mask]
        st.dataframe(filtered, use_container_width=True)
        st.download_button("Download filtered evidence TSV", filtered.to_csv(sep="\t", index=False), "aree_filtered_evidence.tsv")
    else:
        st.info("Run `aree harmonize-demo` to generate evidence.")

with tabs[2]:
    st.subheader("Candidate Biomarker Evidence")
    if scores_path.exists():
        scores = pd.read_csv(scores_path, sep="\t")
        st.dataframe(scores, use_container_width=True)
        selected = st.selectbox("Evidence card", sorted(scores["candidate_id"].astype(str)))
        card = ROOT / "reports" / "evidence_cards" / "{}.md".format(selected.replace("/", "_"))
        if card.exists():
            st.markdown(card.read_text())
    else:
        st.info("Run `aree build-evidence-cards` to generate candidates.")

with tabs[3]:
    st.subheader("Pipeline Status")
    status = pd.DataFrame(
        [
            {"layer": "Study registry", "status": "complete and runnable"},
            {"layer": "Processed harmonization", "status": "complete and runnable for MVP schema"},
            {"layer": "Raw RNA-seq Nextflow", "status": "scaffolded"},
            {"layer": "Raw methylation Nextflow", "status": "scaffolded"},
            {"layer": "Proteomics/metabolomics workflows", "status": "scaffolded"},
            {"layer": "Meta-analysis and scoring", "status": "complete and runnable"},
            {"layer": "Evidence cards", "status": "complete and runnable"},
        ]
    )
    st.dataframe(status, use_container_width=True)

