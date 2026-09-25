# Interpreting Evidence

AREE separates resilience-associated, stress-response, disease-associated, environmental exposure, and suggestive evidence. This distinction prevents overclaiming.

## Meta-analysis

Comparable effects are pooled by feature, feature type, effect-size type, phenotype, and stressor; effects on different scales are never pooled. Effects without a standard error cannot be pooled; they are kept in the meta-analysis table with a `pooling_status` (`pooled`, `single_effect`, or `no_standard_errors`) and counted in `n_effects_excluded`/`excluded_study_ids`, and the report lists per-study coverage. Random-effects summaries report heterogeneity so contradictory findings remain visible.

## Candidate Scores

Candidate scores are prioritization aids. High-scoring candidates should be treated as targets for follow-up validation, not confirmed biomarkers.

## Contradictory Findings

Contradictions are retained in evidence rows, meta-analysis summaries, and evidence cards. A context-dependent candidate can be biologically important but needs careful validation.

