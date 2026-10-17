# Methods

## From an experiment to a gene result

AREE first records what was compared: organism, tissue, life stage, treatment, control, sampling time, replication and any measured performance outcome. Public-study RNA-sequencing data are reanalysed with study-specific designs. Salmon estimates RNA abundance against a documented reference; tximport summarizes transcript estimates at gene level, and DESeq2 estimates expression differences while accounting for the factors specified in each study's model. The scripts and study pages record exclusions, quality checks and exceptions such as processing for 3′ tag sequencing.

Results are put into a common table, a step called **harmonization**. Each row preserves the original gene identifier, reference identifier, comparison, effect estimate, uncertainty, statistical evidence and study limitations. Genome and annotation versions, software methods, source files and file checksums document where the result came from. A shared gene identifier makes records easier to find together; it does not make their experimental comparisons biologically equivalent.

## Reading the results

| Term on the website | Meaning for interpretation |
| --- | --- |
| BioProject | A public collection of sequencing data. It can contain several experimental comparisons. |
| Contrast | The particular treatment-versus-control comparison being analysed. |
| Evidence record | One gene's result in one released comparison. The same gene can have records in several comparisons. A record does not necessarily show a statistically supported difference. |
| Log2 fold change | The estimated expression ratio on a base-2 logarithmic scale. +1 means twice as much RNA in treatment as control; −1 means half as much. Positive does not mean beneficial. |
| Standard error | Uncertainty in the estimated effect under the fitted model, on the same scale as that effect. It depends on replication and model assumptions. |
| Adjusted p-value / FDR | Statistical evidence after accounting for testing many genes. The displayed threshold is adjusted p < 0.05. It is not the probability that a gene predicts survival; FDR control depends on valid model assumptions. Missing values mean no value was available for that test. |
| Biological replication | Independent biological units supporting the comparison. Multiple oysters in one tank or multiple sequencing runs of one RNA sample do not automatically supply independent treatment replicates. Read the study's counted units. |
| RefSeq description | A reference annotation that helps identify possible gene function. It can be predicted or uncharacterized and does not establish that function in this experiment. |
| Evidence card | A summary collecting the available real-study observations for a selected gene, including nonsignificant and exploratory observations. |
| Exploratory / nominal only | A result retained for description whose uncertainty should not be used as reliable inferential evidence. PRJNA735889 currently lacks a model for oysters clustered within tanks. |

## Published real-study findings

The primary website derives findings directly from the committed real harmonized release and its versioned RefSeq annotation tables. It selects **up to ten genes per comparison** for browsing: inference-eligible genes with adjusted p < 0.05, ordered by adjusted p-value, absolute effect, and feature identifier. If fewer than ten qualify, all qualifying genes are displayed. The exploratory PRJNA735889 analysis instead supplies ten descriptive effects ordered by absolute effect, without using nominal significance for selection. Full comparison tables remain downloadable.

Real evidence cards cover the union of these selections and retain every available real-study observation for each feature, including nonsignificant effects. They display study-specific contrasts, uncertainty, phenotype interpretation, biological replication, limitations and input provenance. The selection identifies browsing examples, not the best validated resilience markers.

The primary website does not apply a global biomarker score or pool effects across studies. Differences in exposure, tissue, age, genetic background and shared samples require biological and statistical review. Comparisons withheld from the release remain registered with their reasons; completed reanalysis does not imply an interpretable treatment effect.

## Analysis tools and the synthetic demonstration

The software includes random-effects meta-analysis, which estimates a shared effect while allowing effects to vary between studies, and a transparent candidate score for prioritizing follow-up. These tools are demonstrated with synthetic data in the separate demo section.

The meta-analysis engine groups records by gene or other feature, feature type, effect-size type, phenotype and stressor. It excludes effects without usable standard errors and the flagged tank-dependent exploratory results. It reports variation between effects using Q, tau-squared and I-squared. When multiple usable effects share a replication group, the engine withholds the pooled estimate and heterogeneity as `dependent_effects_require_covariance`. It does not fit covariance between shared samples. Grouping rules also do not establish biological comparability; both issues require review before real synthesis.

Candidate scoring considers independent study groups, reported sample counts, effect magnitude, adjusted significance, direction, phenotype relevance, tissue and life-stage breadth, assay diversity, identifier-mapping confidence, quality flags and heterogeneity. It summarizes magnitude and direction within each effect-size type and recomputes heterogeneity from its inputs. The score is a configurable prioritization aid, not a probability of validation. See [Interpreting evidence](interpreting-evidence.md) for dependence groups and limitations.
