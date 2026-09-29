# AREE Demo Report

This report uses simulated demo evidence to demonstrate the MVP evidence-generation workflow.

## Registered Studies

| study_id | doi | citation | species | assay_type | tissue | life_stage | stressor_class | phenotype | resilience_classification | sample_size | data_status | analysis_status | quality_control_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CGIG_HEAT_RNASEQ_001 | nan | Simulated C. gigas heat challenge RNA-seq study for AREE MVP. | Crassostrea gigas | rnaseq | gill | juvenile | temperature | thermal_tolerance | resilience_associated | 24 | simulated_public_demo | demo_harmonized | simulated_pass |
| CGIG_HEAT_RNASEQ_002 | nan | Simulated C. gigas family heat survival RNA-seq study for AREE MVP. | Crassostrea gigas | rnaseq | mantle | adult | temperature | survival | resilience_associated | 20 | simulated_public_demo | demo_harmonized | simulated_pass |
| CGIG_OA_METHYL_001 | nan | Simulated C. gigas ocean acidification methylation study for AREE MVP. | Crassostrea gigas | methylation | mantle | juvenile | ocean_acidification_pH | acidification_tolerance | resilience_associated | 18 | simulated_public_demo | demo_harmonized | simulated_pass |
| CGIG_PATH_PROTEO_001 | nan | Simulated C. gigas pathogen challenge proteomics study for AREE MVP. | Crassostrea gigas | proteomics | hemolymph | adult | pathogen_challenge | disease_resistance | disease_associated | 16 | simulated_processed_only | demo_harmonized | processed_qc_available |
| CGIG_SALINITY_METAB_001 | nan | Simulated C. gigas salinity stress metabolomics study for AREE MVP. | Crassostrea gigas | metabolomics | whole_body | juvenile | salinity | salinity_tolerance | resilience_associated | 14 | simulated_processed_only | demo_harmonized | simulated_pass_with_annotation_uncertainty |
| CGIG_LARVAL_PROCESSED_001 | nan | Simulated C. gigas larval performance processed-results-only study for AREE MVP. | Crassostrea gigas | processed_results_only | larvae | larva | multi_stressor_exposure | larval_viability | suggestive | 12 | simulated_processed_only | demo_harmonized | limited_metadata |
| CGIG_HEAT_RNASEQ_PRJNA516762 | 10.1186/s12864-019-6003-8 | Liu Y, Li L, Huang B, Wang W, Zhang G. RNAi based transcriptome suggests genes potentially regulated by HSF1 in the Pacific oyster Crassostrea gigas under thermal stress. BMC Genomics. 2019;20:639. | Crassostrea gigas | rnaseq | gill | juvenile | temperature | thermal_tolerance | stress_response | 30 | public_raw_aree_reanalysis_complete | raw_reanalysis_harmonized | All 154,861,045 read pairs quantified with Salmon 2.3.4 against GCF_963853765.1/RS_2024_06; mapping rates 73.85-75.95%; 99 of 22,397 tested genes reach FDR < 0.05 (57 heat-higher, 42 heat-lower); PC1 explains 26.3% of variance and separates most heat from control libraries though heat_A overlaps the controls |
| CGIG_THERMOTOL_RNASEQ_PRJNA694496 | 10.3389/fphys.2021.663023 | Tan Y, Cong R, Qi H, Wang L, Zhang G, Pan Y, Li L. Transcriptomics Analysis and Re-sequencing Reveal the Mechanism Underlying the Thermotolerance of an Artificial Selection Population of the Pacific Oyster. Front Physiol. 2021;12:663023. | Crassostrea gigas | rnaseq | whole_larva | umbo_larva | temperature | thermal_tolerance | resilience_associated | 6 | public_raw_aree_reanalysis_complete | raw_reanalysis_harmonized | All 132,413,989 read pairs quantified with Salmon 1.10.3; mapping rates 42.13-63.94%; PC1 explains 90.6% of variance and separates selected from control; selected_A is displaced on PC2 and has six genes above the Cook's-distance threshold, but leave-one-out effects correlate 0.962 with the primary model and preserve the direction of all 8,802 primary significant genes |
| CGIG_OSHV1_RNASEQ_PRJNA1329250 | 10.1016/j.fsi.2026.111154 | Calla B, Thompson NF, Burge CA. Population-specific transcriptomics of Pacific oyster after exposure to a highly pathogenic, globally distributed virus. Fish Shellfish Immunol. 2026;171:111154. | Crassostrea gigas | rnaseq | whole_soft_tissue | spat | pathogen_challenge | immune_responsiveness | disease_associated | 42 | public_raw_aree_reanalysis_complete | raw_reanalysis_harmonized | All 1,564,162,034 read pairs quantified with Salmon 2.3.4 against GCF_963853765.1/RS_2024_06; mapping rates 44.24-61.66% (median 55.93%), well below the 73.9-76.0% of PRJNA516762 but uniform across arms (control 57.0% versus challenged 55.4%, Welch p = 0.30) and across both lineages, so it reflects library preparation rather than viral reads displacing host reads; 5,446 of 25,734 tested genes reach FDR < 0.05 (3,237 challenge-higher, 2,209 challenge-lower); PC1 explains 59.1% of variance and orders challenged above control though six control libraries fall inside the challenged range, and PC2 (16.6%) tracks the Midori/Miyagi split that ~ population + condition blocks on |
| CGIG_OA_RNASEQ_PRJNA826964 | nan | Institute of Oceanology, Chinese Academy of Sciences. The compromised energy management of Pacific oysters (Crassostrea gigas) under ocean acidification conditions. NCBI BioProject PRJNA826964; public 2022. No linked publication identified. | Crassostrea gigas | rnaseq | hepatopancreas | not_recorded | ocean_acidification_pH | metabolic_resilience | stress_response | 6 | public_raw_aree_reanalysis_complete | raw_reanalysis_complete_not_harmonized | All 435,657,730 read pairs quantified with Salmon 2.3.4 against GCF_963853765.1/RS_2024_06; mapping rates 78.62-81.26% (median 79.79%), equal across arms (control 79.9% versus acidified 79.9%); expressed-SNP allele frequencies at 2,837 sites show the three libraries in each day-by-condition cell are one pooled RNA sample (within-cell mean difference 0.040-0.057 against 0.039 expected from resequencing the same RNA, r = 0.95; between cells 0.078-0.317), so the 18 libraries were collapsed to 6 pools and tested as ~ timepoint + condition; 890 of 24,681 tested genes reach FDR < 0.05 (301 acidified-higher, 589 acidified-lower) against 10,797 when the libraries were treated as replicates; PC1 (35.9%) separates acidified_d56 from the other five pools and PC2 (29.7%) separates the arms |
| CGIG_SALINITY_RNASEQ_PRJNA756710 | 10.3389/fimmu.2022.859975 | Li X, Yang B, Shi C, Wang H, Yu R, Li Q, Liu S. Synergistic interaction of low salinity stress with Vibrio infection causes mass mortalities in the oyster by inducing host microflora imbalance and immune dysregulation. Front Immunol. 2022;13:859975. | Crassostrea gigas | rnaseq | gill | not_recorded | salinity | disease_susceptibility | disease_associated | 12 | public_raw_aree_reanalysis_complete | raw_reanalysis_harmonized | All 260,215,874 read pairs quantified with Salmon 2.3.4 against GCF_963853765.1/RS_2024_06; mapping rates 70.78-80.24% (median 76.97%), consistent with the 73.9-82.5% the publication reports against GCA_902806645.1 and flat across arms (30 ppt 77.7% versus 10 ppt 76.6%, Welch p = 0.49); 3,045 of 21,098 tested genes reach FDR < 0.05 (1,825 low-salinity-higher, 1,220 lower); PC1 explains 27.9% of variance and separates all six 10 ppt libraries from all six 30 ppt libraries; low_salinity_h12_2 is displaced on PC2 and has the most genes above the Cook's-distance threshold (84), but leaving it out gives effects correlated 0.929 with the primary model, preserves the direction of all 3,045 primary significant genes and retains 2,618 of them at FDR < 0.05 |
| CGIG_STARVATION_RNASEQ_PRJEB86618 | 10.1016/j.isci.2025.114333 | Duperret L, Valdivieso A, Petton B, Morga B, Faury N, de Lorgeril J, Corporeau C, Degremont L, Allienne JF, Henry S, Turtoi A, Pouzadoux J, Romatif O, Courtay G, Monaco CJ, Toulza E, Lagorce A, Pernet F, Vidal-Dupiol J, Mitta G. Food deprivation enhances disease resistance: Underlying mechanisms in oysters confronted with pacific oyster mortality syndrome. iScience. 2026;29(1):114333. | Crassostrea gigas | rnaseq | whole_soft_tissue | spat | nutritional_limitation | disease_resistance | resilience_associated | 12 | public_raw_aree_reanalysis_complete | raw_reanalysis_harmonized | All 318,481,607 single-end reads quantified with Salmon 2.3.4 (sequence-bias correction; GC-bias correction is not fitted for single-end reads) against GCF_963853765.1/RS_2024_06; mapping rates 65.35-69.05% (median 67.26%), below the publication's 82.3% STAR genome alignment as expected for transcriptome-only quantification, and about two points lower in starved than fed pools (66.2% versus 68.2%, Welch p = 0.0003, with decoy-genome rates 9.8% versus 8.0%); 11,930 of 25,253 tested genes reach FDR < 0.05 (5,941 starved-higher, 5,989 starved-lower); PC1 (49.1%) and PC2 (39.8%) separate the four diet-by-family cells; within-cell BCV is 9.8%, below the 17-38% of the other AREE oyster studies, which fits pools of 10 full sibs. The pools are independent: expressed-SNP allele frequencies at 477 sites differ within cells by 0.098 against 0.042 simulated for one RNA resequenced, and no library shows single-animal genotype peaks |

## Evidence Counts

| phenotype | stressor | feature_type | records |
| --- | --- | --- | --- |
| acidification_tolerance | ocean_acidification_pH | genomic_region | 2 |
| disease_resistance | pathogen_challenge | protein | 2 |
| larval_viability | multi_stressor_exposure | gene | 2 |
| salinity_tolerance | salinity | metabolite | 2 |
| survival | temperature | gene | 3 |
| thermal_tolerance | temperature | gene | 4 |

## Meta-analysis Coverage

Only effects with a standard error can be pooled; effects without one are listed here rather than dropped silently.

| study_id | effects | poolable_effects | not_poolable |
| --- | --- | --- | --- |
| CGIG_HEAT_RNASEQ_001 | 4 | 4 | 0 |
| CGIG_HEAT_RNASEQ_002 | 3 | 3 | 0 |
| CGIG_LARVAL_PROCESSED_001 | 2 | 2 | 0 |
| CGIG_OA_METHYL_001 | 2 | 2 | 0 |
| CGIG_PATH_PROTEO_001 | 2 | 2 | 0 |
| CGIG_SALINITY_METAB_001 | 2 | 2 | 0 |

## Candidate Scores

| candidate_id | score | category | n_studies | total_biological_sample_size | assay_diversity | direction_consistency | consistency_flag | best_adjusted_p_value | mean_mapping_confidence_score | known_limitations |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:LOC105317001 | 0.6576 | High-priority cross-study candidate | 3 | 56 | 1 | 1.0 | reasonably consistent | 0.008 | 1.0 | none; processed_only_limited_metadata |
| NCBI:LOC105317002 | 0.6128 | Multi-omics convergence candidate | 3 | 60 | 2 | 0.5 | conflicting or context-dependent | 0.018 | 0.933 | conflicting_direction; none |
| NCBI:LOC105317004 | 0.5807 | Multi-omics convergence candidate | 3 | 54 | 3 | nan | not replicated within an effect-size type | 0.07 | 0.467 | none; processed_only |
| NCBI:LOC105317003 | 0.5803 | Multi-omics convergence candidate | 3 | 54 | 2 | nan | not replicated within an effect-size type | 0.03 | 0.8 | nearest_gene_annotation; none; suggestive_phenotype |
| CGI_99999 | 0.4073 | Emerging candidate requiring replication | 1 | 24 | 1 | nan | not replicated within an effect-size type | 0.09 | 0.2 | imperfect_identifier_mapping |
| glutathione_related_feature | 0.3832 | Emerging candidate requiring replication | 1 | 14 | 1 | nan | not replicated within an effect-size type | 0.04 | 0.2 | annotation_inferred |
| unknown_metabolite_404 | 0.3734 | Emerging candidate requiring replication | 1 | 14 | 1 | nan | not replicated within an effect-size type | 0.14 | 0.0 | unresolved_metabolite |

## Interpretation Guardrail

Candidates are prioritized associations. A statistically significant single-study result is not treated as validation.
