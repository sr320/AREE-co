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
| CGIG_HEAT_RNASEQ_PRJNA516762 | 10.1186/s12864-019-6003-8 | Liu Y, Li L, Huang B, Wang W, Zhang G. RNAi based transcriptome suggests genes potentially regulated by HSF1 in the Pacific oyster Crassostrea gigas under thermal stress. BMC Genomics. 2019;20:639. | Crassostrea gigas | rnaseq | gill | juvenile | temperature | thermal_tolerance | stress_response | 30 | public_raw_and_processed | processed_results_harmonized | Publication-reported QC and 75.18-76.51% alignment; not independently reanalyzed |
| CGIG_THERMOTOL_RNASEQ_PRJNA694496 | 10.3389/fphys.2021.663023 | Tan Y, Cong R, Qi H, Wang L, Zhang G, Pan Y, Li L. Transcriptomics Analysis and Re-sequencing Reveal the Mechanism Underlying the Thermotolerance of an Artificial Selection Population of the Pacific Oyster. Front Physiol. 2021;12:663023. | Crassostrea gigas | rnaseq | whole_larva | umbo_larva | temperature | thermal_tolerance | resilience_associated | 6 | public_raw_aree_reanalysis_complete | raw_reanalysis_harmonized | All 132,413,989 read pairs quantified with Salmon 1.10.3; mapping rates 42.13-63.94%; PC1 explains 90.6% of variance and separates selected from control; selected_A is displaced on PC2 and has six genes above the Cook's-distance threshold, but leave-one-out effects correlate 0.962 with the primary model and preserve the direction of all 8,802 primary significant genes |

## Evidence Counts

| phenotype | stressor | feature_type | records |
| --- | --- | --- | --- |
| acidification_tolerance | ocean_acidification_pH | genomic_region | 2 |
| disease_resistance | pathogen_challenge | protein | 2 |
| larval_viability | multi_stressor_exposure | gene | 2 |
| salinity_tolerance | salinity | metabolite | 2 |
| survival | temperature | gene | 3 |
| thermal_tolerance | temperature | gene | 4 |

## Candidate Scores

| candidate_id | score | category | n_studies | total_biological_sample_size | assay_diversity | direction_consistency | consistency_flag | best_adjusted_p_value | mean_mapping_confidence_score | known_limitations |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:LOC105317001 | 0.6624 | High-priority cross-study candidate | 3 | 56 | 1 | 1.0 | reasonably consistent | 0.008 | 1.0 | none; processed_only_limited_metadata |
| NCBI:LOC105317002 | 0.6333 | Multi-omics convergence candidate | 3 | 60 | 2 | 0.667 | conflicting or context-dependent | 0.018 | 0.933 | conflicting_direction; none |
| NCBI:LOC105317004 | 0.6029 | Multi-omics convergence candidate | 3 | 54 | 3 | 0.667 | conflicting or context-dependent | 0.07 | 0.467 | none; processed_only |
| NCBI:LOC105317003 | 0.6026 | Multi-omics convergence candidate | 3 | 54 | 2 | 0.667 | conflicting or context-dependent | 0.03 | 0.8 | nearest_gene_annotation; none; suggestive_phenotype |
| CGI_99999 | 0.474 | Emerging candidate requiring replication | 1 | 24 | 1 | 1.0 | reasonably consistent | 0.09 | 0.2 | imperfect_identifier_mapping |
| glutathione_related_feature | 0.4498 | Emerging candidate requiring replication | 1 | 14 | 1 | 1.0 | reasonably consistent | 0.04 | 0.2 | annotation_inferred |
| unknown_metabolite_404 | 0.44 | Emerging candidate requiring replication | 1 | 14 | 1 | 1.0 | reasonably consistent | 0.14 | 0.0 | unresolved_metabolite |

## Interpretation Guardrail

Candidates are prioritized associations. A statistically significant single-study result is not treated as validation.
