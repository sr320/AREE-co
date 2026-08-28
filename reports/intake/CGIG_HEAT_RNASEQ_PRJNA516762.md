# Real-study intake: PRJNA516762

## Scope

- Study: Liu et al. (2019), *RNAi based transcriptome suggests genes potentially regulated by HSF1 in the Pacific oyster Crassostrea gigas under thermal stress*.
- DOI: <https://doi.org/10.1186/s12864-019-6003-8>
- Raw-data accessions: BioProject `PRJNA516762`; SRA study `SRP181984`; 12 paired-end Illumina HiSeq X Ten runs.
- AREE contrast: untreated control versus 35 C heat shock for 2 h, excluding the HSF1 RNAi arms.
- Biological replication for this contrast: 3 pooled libraries per condition, each pooling 5 oysters (30 oysters represented in 6 libraries).

## Processed-result QC

- The publication reports 150 significant genes for the direct heat-versus-control contrast; the source workbook and derived AREE table both contain 150 records.
- Direction counts: 76 upregulated and 74 downregulated.
- All reported directions agree with the sign of the published log2 fold change.
- FDR range: `3.0261e-86` to `0.049902`; all records satisfy the publication's FDR threshold.
- Raw p-values and standard errors are not available in the supplement and remain null in AREE.
- All 150 legacy `CGI_*` or newly assembled identifiers remain unresolved against the current demo mapping and must not yet be treated as cross-study standardized features.

## Largest published effects

| Direction | Feature | log2 fold change | FDR |
| --- | --- | ---: | ---: |
| Up | `CGI_10002387` | 4.3219 | 3.0261e-86 |
| Up | `CGI_10002823` | 2.0048 | 1.1331e-10 |
| Up | `CGI_10009093` | 1.8845 | 2.3745e-14 |
| Up | `CGI_10002068` | 1.8766 | 1.5230e-09 |
| Up | `CGI_10002950` | 1.8217 | 1.2167e-14 |
| Down | `CGI_10020949` | -3.5236 | 1.9524e-37 |
| Down | `CGI_10018182` | -2.7578 | 7.1855e-24 |
| Down | `CGI_10001445` | -2.5893 | 1.8609e-18 |
| Down | `CGI_10010392` | -2.4298 | 1.0999e-22 |
| Down | `CGI_10016700` | -2.3808 | 1.1136e-15 |

## Interpretation guardrail

These are significant molecular heat-response associations from a pooled, single-study contrast. The experiment did not directly measure survival or organism-level thermal tolerance, so these records are classified as `stress_response`, not validated resilience biomarkers.

## Next analytical step

Build a versioned mapping from the legacy `CGI_*` identifiers used with assembly `GCA_000297895.1` to the current *C. gigas* gene annotation. After mapping, curate an independent thermal-challenge study with compatible effect estimates. Random-effects meta-analysis should wait until at least two independent studies have standardized identifiers and usable uncertainty estimates.
