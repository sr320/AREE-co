# Repository and website review — 17 October 2026

Reviewed the repository at `1b0b2f3` against the published [AREE website](https://sr320.github.io/AREE-co/), then revised its documentation and website-generation text. The initially clean local checkout was 31 commits behind the already-fetched `origin/main`; it was fast-forwarded to the revision identified on the live page before editing.

## Accuracy checks

The live homepage and all 31 rows of its downloadable study summary agree with the rebuilt site and committed evidence:

| Item | Verified count |
| --- | ---: |
| Public BioProjects represented by registered comparisons | 20 |
| Registered real-study comparisons | 31 |
| Comparisons with released harmonized evidence | 24 |
| Gene-by-comparison evidence records | 518,407 |
| Displayed gene-by-comparison findings | 238 |
| Distinct real gene evidence cards | 200 |

Every study YAML file passes the study schema. The site loader checks registry CSV/YAML consistency, final sample counts, study identifiers and duplicate genes. The review additionally verified all 24 processed-input SHA-256 checksums and compared every released gene's effect, standard error, p-value and adjusted p-value against its processed input, without discrepancies. These checks establish internal consistency and file provenance, rather than independently re-estimating the underlying biological effects.

## Corrections made

- Replaced outdated introductory descriptions that framed the project as a simulated-data demonstration awaiting public-study curation.
- Corrected “ten findings per study” to “up to ten”: the tire-particle comparison currently has eight qualifying genes. The 238 displayed findings therefore need not equal ten times 24.
- Clarified that the PRJNA694496 report's historical overlap uses a publication-derived PRJNA516762 gene list. PRJNA516762 now also has a later raw-data evidence release. Updated the report generator to describe its supplied input without always calling it a published significant-gene list.
- Corrected the architecture description to distinguish study-specific scripts used for public reanalyses from the reusable paired-end Nextflow pipeline.
- Added an accessible explanation of the work to the homepage and README, expanded the biological overview, and added a methods glossary linked from findings and evidence cards.
- Explained what project, comparison, evidence-record and card counts represent; how to read fold changes and uncertainty; and why RNA differences do not establish improved survival or causality.
- Kept the real-study/synthetic-demo distinction and the existing restrictions around shared controls, exploratory tank-dependent inference and cross-study pooling explicit.

## Verification and limits

The Python suite completed with 163 tests passing and four local-HTTP fixture setup errors caused by sandbox socket restrictions. Rerunning the affected download-helper module outside the sandbox passed all nine tests, including the four affected cases. Ruff and `git diff --check` pass. Generated Markdown/Quarto pages have no missing relative link targets.

The complete Quarto site renders successfully: 300 HTML pages, with no broken local page links or fragment anchors. These changes do not rerun public raw-read analyses, independently validate every publication's biological claims, or validate resilience biomarkers. Existing study-design limitations remain, including the need for a tank-aware PRJNA735889 analysis and careful interpretation of genetic background, developmental stage and shared samples. The reanalysis-complete count includes comparisons withheld from the released evidence.

Edits are local source changes. The live website will change only after they are committed, pushed to `main` and successfully deployed by the Pages workflow.
