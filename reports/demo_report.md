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
| NCBI:LOC105317001 | 0.6955 | High-priority cross-study candidate | 3 | 56 | 1 | 1.0 | reasonably consistent | 0.008 | 1.0 | none; processed_only_limited_metadata |
| NCBI:LOC105317002 | 0.665 | Multi-omics convergence candidate | 3 | 60 | 2 | 0.667 | conflicting or context-dependent | 0.018000000000000002 | 0.9329999999999999 | conflicting_direction; none |
| NCBI:LOC105317004 | 0.633 | Multi-omics convergence candidate | 3 | 54 | 3 | 0.667 | conflicting or context-dependent | 0.07 | 0.467 | none; processed_only |
| NCBI:LOC105317003 | 0.6327 | Multi-omics convergence candidate | 3 | 54 | 2 | 0.667 | conflicting or context-dependent | 0.03 | 0.8 | nearest_gene_annotation; none; suggestive_phenotype |
| CGI_99999 | 0.4977 | Emerging candidate requiring replication | 1 | 24 | 1 | 1.0 | reasonably consistent | 0.09 | 0.2 | imperfect_identifier_mapping |
| glutathione_related_feature | 0.4723 | Emerging candidate requiring replication | 1 | 14 | 1 | 1.0 | reasonably consistent | 0.04 | 0.2 | annotation_inferred |
| unknown_metabolite_404 | 0.462 | Emerging candidate requiring replication | 1 | 14 | 1 | 1.0 | reasonably consistent | 0.14 | 0.0 | unresolved_metabolite |

## Interpretation Guardrail

Candidates are prioritized associations. A statistically significant single-study result is not treated as validation.
