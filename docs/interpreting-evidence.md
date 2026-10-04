# Interpreting Evidence

AREE separates resilience-associated, stress-response, disease-associated, environmental exposure, and suggestive evidence. This distinction prevents overclaiming.

## How expression findings relate to resilience

**Differentially expressed genes are potential clues to resilience; differential expression alone does not establish a resilience marker.** It shows that RNA abundance differs between the groups being compared. The biological interpretation depends on the comparison and whether resilience was measured.

In aquaculture, resilience needs an observable outcome, such as survival during a challenge, maintained growth, lower disease burden, or recovery afterward. Expression becomes evidence for resilience when it is associated with that outcome. A marker can track an outcome without causing it; a causal claim requires additional experiments. The [FDA–NIH biomarker framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC5813875/) distinguishes biological indicators from biomarkers associated with outcomes and specific uses.

| Comparison | What the expression difference supports |
| --- | --- |
| Heat-exposed versus unexposed oysters | A **heat-response marker**. The difference could reflect protection, damage, or both. |
| A thermotolerance-selected population versus controls, with better survival demonstrated | A **candidate resilience-associated marker**. Other differences between the populations could also explain expression. |
| Expression measured before challenge that consistently forecasts subsequent survival in independent animals | Stronger evidence for a **marker useful for forecasting resilience** in that tested context. This is a validation target, not a claim about the current cards. |

The current studies illustrate these distinctions. PRJNA516762 compares heat exposure with an untreated control; its released contrast measures a molecular heat response without an organism-level thermal-tolerance outcome. PRJNA694496 compares larvae from a thermotolerance-selected population with controls at ambient conditions. The selected population had better heat-challenge and summer survival in the [original study](https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2021.663023/full). AREE interprets its expression differences as candidates associated with a more tolerant population; they have not been shown to predict each individual oyster's survival.

## Expression direction and statistical significance

**Higher expression does not automatically mean greater resilience.** For example, high expression of a heat-shock gene could indicate an effective protective response or severe cellular stress. A tolerant oyster might need less induction because it experiences less damage. This is an illustrative possibility, not a conclusion about a particular gene in the evidence cards.

Sampling time changes the interpretation: baseline preparedness, acute response, and recovery represent different biological states. Tissue, life stage, family, exposure intensity and cell composition can also affect measured expression. On the findings pages, a positive effect means higher expression in the stated treatment relative to its control; it does not mean a favorable effect on resilience.

A small adjusted p-value supports an expression difference under the statistical model. It does not measure how well a gene forecasts survival or establish that changing its expression improves tolerance. Large expression changes and plausible gene functions can guide follow-up, but they do not replace outcome evidence.

## From a finding to a useful resilience marker

1. **Connect expression to a measured outcome**, ideally in the same animals or clearly defined families and experimental units. Specify whether the intended marker measures current response, future performance, or recovery.
2. **Validate the relationship independently**, with tissue, life stage, exposure and sampling time specified and experimental dependence accounted for.
3. **Test practical performance**, including whether expression adds useful information about outcomes beyond factors such as family, size and environment. An assay intended for use before a challenge must be validated at that sampling time.
4. **Test causality separately** if claiming that a gene contributes to resilience, rather than serving as an indicator.

The current real evidence cards are **study-specific molecular findings with varying relevance to resilience**. Read their phenotype interpretation and study limitations to see whether survival, performance or recovery was measured and how it relates to the expression samples. The ten-per-study selection is a browsing choice based on statistical evidence or descriptive effect size; it does not identify the ten best resilience markers. Presence in multiple cards or studies is not sufficient validation when contrasts differ or samples are shared.

## Meta-analysis

The primary website presents real study-level findings without cross-study pooling. The following describes the analysis engine; synthetic pooled results are available in the secondary demo section.

Comparable effects are pooled by feature, feature type, effect-size type, phenotype, and stressor; effects on different scales are never pooled. Effects without usable uncertainty or marked as exploratory because tank clustering is unmodeled cannot be pooled. They remain in the meta-analysis table with a `pooling_status` (`pooled`, `single_effect`, `no_standard_errors`, or `no_inference_eligible_effects`) and are counted in `n_effects_excluded`/`excluded_study_ids`; the report lists per-study coverage. Random-effects summaries report heterogeneity so contradictory findings remain visible.

## Candidate Scores

The primary real-study findings and cards do not use a global candidate score. Synthetic scores demonstrate the scoring workflow in the secondary demo section.

Candidate scores are prioritization aids. High-scoring candidates should be treated as targets for follow-up validation, not confirmed biomarkers.

`n_studies` counts registered study IDs; `n_independent_studies` counts replication groups. Evidence with a shared `dependence_group=GROUP_ID` quality flag counts as one group for the replication reward and ranking category. Directional consistency is computed after averaging correlated contrasts within each group and effect-size type. For shared samples, `total_biological_sample_size` is a conservative lower bound (the largest study total per group), identified by `sample_size_status`; it is not the exact number of unique animals or libraries. Reported sample sizes may describe animals or pooled libraries, so biological replication must also be checked on the study detail page.

PRJEB86646 and PRJEB86618 share six DECICOMP control libraries and carry the same dependence group. Their agreement alone cannot establish independent cross-study replication.

The two PRJNA913164 records (30 C seawater, and 30 C followed by 44 C emersion) share the same 23 diploid and triploid 20 C control oysters and carry the dependence group `PRJNA913164_20C_controls`, so they too count as one replication group.

PRJNA735889 currently uses an oyster-level model on 50 oysters from ten tanks without accounting for tank clustering. Its p-values, FDR values and standard errors are exploratory until a tank-aware reanalysis is available. Evidence marked `exploratory_tank_clustering_unmodeled` receives no significance reward, cannot promote a candidate to the high-priority cross-study category, and remains visible with its limitations. `best_adjusted_p_value` in scores excludes these exploratory values; original nominal values remain in the evidence for provenance.

## Contradictory Findings

Contradictions are retained in evidence rows, meta-analysis summaries, and evidence cards. A context-dependent candidate can be biologically important but needs careful validation.
