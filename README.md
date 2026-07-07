# Aquaculture Resilience Evidence Engine (AREE)

AREE is an open, reproducible evidence-generation system for aquaculture resilience biomarker discovery. The MVP focuses on Pacific oyster (*Crassostrea gigas*) and is structured so additional shellfish species and aquaculture organisms can be added without changing the core model.

AREE is not a static list of papers. It provides a versioned path from public study registration, through harmonized processed-result intake or raw workflow scaffolds, into cross-study evidence tables, meta-analysis summaries, transparent candidate prioritization, evidence cards, and user-facing reports.

## What Is Runnable Now

- Study registration validation against JSON Schema.
- Registry ingestion from YAML into `registry/study_registry.csv`.
- Processed-result harmonization for RNA-seq, methylation, proteomics, and metabolomics-style tables.
- Demo cross-study random-effects meta-analysis.
- Transparent biomarker candidate scoring.
- Markdown evidence-card generation.
- Demo report generation.
- A Streamlit browser for studies, evidence, and candidate cards.
- Nextflow workflow scaffolds for raw-data reanalysis.

Synthetic demo data are clearly labeled as simulated and are intended to exercise the full system shape before raw public datasets are curated.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,app]"
```

If you prefer to use an existing Python environment:

```bash
python3 -m pip install -e ".[dev,app]"
```

## Demo Commands

```bash
aree validate-study registry/studies/CGIG_HEAT_RNASEQ_001.yaml
aree register-study registry/studies/CGIG_HEAT_RNASEQ_001.yaml
aree harmonize --study CGIG_HEAT_RNASEQ_001 --input data/demo/processed/CGIG_HEAT_RNASEQ_001_rnaseq.tsv
aree harmonize-demo
aree meta-analyze --phenotype thermal_tolerance --feature-type gene
aree build-evidence-cards --phenotype survival
aree build-demo-report
```

Optional interfaces:

```bash
streamlit run app/main.py
quarto render docs/
```

`quarto render docs/` is scaffolded for users with Quarto installed. The runnable Python report generator is `aree build-demo-report`.

## Repository Tree

```text
AREE/
├── README.md
├── LICENSE
├── CITATION.cff
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
├── docs/
├── schemas/
├── registry/
│   ├── studies/
│   ├── controlled_vocabularies/
│   └── study_registry.csv
├── workflows/
│   ├── rnaseq/
│   ├── methylation/
│   ├── proteomics/
│   └── metabolomics/
├── modules/
├── containers/
├── config/
├── data/
│   ├── demo/
│   ├── reference/
│   └── mappings/
├── src/
│   └── aree/
├── app/
├── reports/
├── tests/
├── notebooks/
└── .github/
    └── workflows/
```

## Objective 1 Component Map

| AREE component | Objective 1 function |
|---|---|
| Study registry | Standardizes public dataset intake and study characterization. |
| Controlled vocabularies | Normalizes phenotype, stressor, assay, tissue, and quality labels. |
| JSON schemas | Makes dataset metadata and evidence tables machine-validatable. |
| Nextflow scaffolds | Defines reproducible raw-data reanalysis entry points. |
| Processed-result harmonizers | Allows studies without raw data to contribute transparent evidence. |
| Identifier mapping | Preserves original IDs while assigning comparable reference identifiers and mapping confidence. |
| Evidence schema | Converts assay-specific outputs into shared biomarker evidence. |
| Meta-analysis code | Pools comparable effects and reports heterogeneity. |
| Candidate scoring | Prioritizes repeated, interpretable, quality-aware biomarker associations without claiming validation. |
| Evidence cards | Produces database-ready summaries with limitations and recommended validation steps. |
| Streamlit/report outputs | Gives researchers and breeders a practical way to inspect and download evidence. |

## Production Status

| Area | Status |
|---|---|
| Study schema validation | Complete and runnable. |
| Demo registry and synthetic processed results | Complete and runnable. |
| Processed-result harmonization | Complete for MVP input schemas; production importers need assay-specific expansion. |
| Meta-analysis | Complete MVP random-effects implementation; should be cross-checked with R `metafor` for production. |
| Candidate scoring | Complete transparent MVP; weights should be reviewed by domain experts. |
| Evidence cards | Complete Markdown output; plots are textual summaries in MVP. |
| Streamlit app | Complete lightweight browser; production search and visualizations remain future work. |
| Raw-data Nextflow workflows | Scaffolded, not production-ready. |
| Containers | Scaffolded. |
| CI | Included; depends on runner availability for optional Quarto/Streamlit dependencies. |

## First Curation Targets

Suggested first public *C. gigas* dataset classes to curate:

1. RNA-seq thermal challenge study with survival or heat-tolerance phenotype.
2. RNA-seq ocean acidification exposure study with growth or shell-formation phenotype.
3. RNA-seq pathogen challenge study with disease resistance or pathogen-load phenotype.
4. WGBS or EM-seq temperature or acidification study with methylation changes linked to resilience.
5. Proteomics heat or pathogen challenge study with tissue-specific abundance changes.
6. Metabolomics hypoxia or salinity exposure study with survival or recovery outcomes.
7. Larval viability study under acidification or temperature stress.
8. Family or breeding-line comparison under thermal challenge.
9. Salinity or freshwater exposure study with recovery following stress.
10. Multi-stressor exposure study combining temperature with pH, hypoxia, or pathogen challenge.

Do not invent accession numbers during curation. Add only accessions verified from primary repositories or publications.

