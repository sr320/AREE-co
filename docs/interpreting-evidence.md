# Interpreting Evidence

AREE separates resilience-associated, stress-response, disease-associated, environmental exposure, and suggestive evidence. This distinction prevents overclaiming.

## Meta-analysis

Comparable effects are pooled by feature, feature type, effect-size type, phenotype, and stressor; effects on different scales are never pooled. Effects without usable uncertainty or marked as exploratory because tank clustering is unmodeled cannot be pooled. They remain in the meta-analysis table with a `pooling_status` (`pooled`, `single_effect`, `no_standard_errors`, or `no_inference_eligible_effects`) and are counted in `n_effects_excluded`/`excluded_study_ids`; the report lists per-study coverage. Random-effects summaries report heterogeneity so contradictory findings remain visible.

## Candidate Scores

Candidate scores are prioritization aids. High-scoring candidates should be treated as targets for follow-up validation, not confirmed biomarkers.

`n_studies` counts registered study IDs; `n_independent_studies` counts replication groups. Evidence with a shared `dependence_group=GROUP_ID` quality flag counts as one group for the replication reward and ranking category. Directional consistency is computed after averaging correlated contrasts within each group and effect-size type. For shared samples, `total_biological_sample_size` is a conservative lower bound (the largest study total per group), identified by `sample_size_status`; it is not the exact number of unique animals or libraries. Reported sample sizes may describe animals or pooled libraries, so biological replication must also be checked on the study detail page.

PRJEB86646 and PRJEB86618 share six DECICOMP control libraries and carry the same dependence group. Their agreement alone cannot establish independent cross-study replication.

PRJNA735889 currently uses an oyster-level model on 50 oysters from ten tanks without accounting for tank clustering. Its p-values, FDR values and standard errors are exploratory until a tank-aware reanalysis is available. Evidence marked `exploratory_tank_clustering_unmodeled` receives no significance reward, cannot promote a candidate to the high-priority cross-study category, and remains visible with its limitations. `best_adjusted_p_value` in scores excludes these exploratory values; original nominal values remain in the evidence for provenance.

## Contradictory Findings

Contradictions are retained in evidence rows, meta-analysis summaries, and evidence cards. A context-dependent candidate can be biologically important but needs careful validation.
