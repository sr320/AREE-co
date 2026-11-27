# Candidate Scoring Weight & Quality-Flag Audit

Evaluated `score_table(..., method='reml', shared_control_rho=0.5)` across all 5 comparability-reviewed strata under equal weights, replication-heavy, outcome/quality-heavy, and leave-one-component-out perturbations.

## Weight Stability Summary

| scheme | mean_spearman | min_spearman | mean_top25_jaccard | mean_top100_jaccard |
| --- | --- | --- | --- | --- |
| equal_weights | 0.9929 | 0.9915 | 0.9264 | 0.9539 |
| replication_heavy_2x | 0.9863 | 0.9818 | 0.9141 | 0.8872 |
| outcome_quality_heavy_2x | 0.992 | 0.9882 | 0.9264 | 0.9349 |
| drop_n_studies | 0.9588 | 0.9189 | 0.8277 | 0.782 |
| drop_total_sample_size | 0.9905 | 0.9615 | 0.9294 | 0.9352 |
| drop_effect_magnitude | 0.9621 | 0.955 | 0.2978 | 0.3777 |
| drop_significance | 0.974 | 0.9511 | 0.3086 | 0.4142 |
| drop_direction_consistency | 0.937 | 0.9129 | 0.7099 | 0.6633 |
| drop_phenotype_relevance | 0.9967 | 0.9953 | 0.918 | 0.9391 |
| drop_context_breadth | 0.9895 | 0.9835 | 0.9294 | 0.9497 |
| drop_assay_diversity | 0.9967 | 0.9953 | 0.918 | 0.9391 |
| drop_mapping_confidence | 0.9991 | 0.9988 | 0.9418 | 0.965 |
| drop_data_quality | 0.9933 | 0.9745 | 0.9294 | 0.9386 |

## Quality-Flag Classification Across All Harmonized Real Studies

| quality_flag | classified_as |
| --- | --- |
| acute_heat_response | provenance_non_penalty |
| constitutive_expression | provenance_non_penalty |
| dependence_group=DECICOMP_PRJEB86525_0h_controls | provenance_non_penalty |
| dependence_group=PRJEB58194_F1_TT | provenance_non_penalty |
| dependence_group=PRJNA298285_ambient_larvae | provenance_non_penalty |
| dependence_group=PRJNA407831_0h_pools | provenance_non_penalty |
| dependence_group=PRJNA515169_mock_pools | provenance_non_penalty |
| dependence_group=PRJNA593309_21C_tanks | provenance_non_penalty |
| dependence_group=PRJNA852380_controls | provenance_non_penalty |
| dependence_group=PRJNA856813_gill_controls | provenance_non_penalty |
| dependence_group=PRJNA913164_20C_controls | provenance_non_penalty |
| deseq2_significance_unavailable | caveat_penalty |
| direct_exposure | provenance_non_penalty |
| duplicate_animal_libraries_collapsed | caveat_penalty |
| during_oshv1_infection | provenance_non_penalty |
| experiment_blocked_contrast | provenance_non_penalty |
| family_blocked_contrast | provenance_non_penalty |
| fdr_lt_0.05 | provenance_non_penalty |
| gametogenesis_blocked_contrast | provenance_non_penalty |
| harmful_algal_bloom_response | provenance_non_penalty |
| heat_24h | provenance_non_penalty |
| heat_6h | provenance_non_penalty |
| infection_background_all_arms | provenance_non_penalty |
| larval_development | provenance_non_penalty |
| lineage_balanced_contrast | provenance_non_penalty |
| low_salinity_response | provenance_non_penalty |
| multiple_stressor_response | provenance_non_penalty |
| no_ambient_ph_reference | caveat_penalty |
| no_linked_publication | caveat_penalty |
| nutritional_limitation_response | provenance_non_penalty |
| ocean_acidification_response | provenance_non_penalty |
| one_oyster_per_tank | provenance_non_penalty |
| one_tank_per_ph_level | caveat_penalty |
| pathogen_challenge_response | provenance_non_penalty |
| pesticide_mixture | provenance_non_penalty |
| ph_band_contrast_tipping_window_excluded | caveat_penalty |
| ploidy_blocked_contrast | provenance_non_penalty |
| pollutant_exposure_response | provenance_non_penalty |
| pooled_libraries | caveat_penalty |
| pooled_rna_samples | caveat_penalty |
| population_blocked_contrast | provenance_non_penalty |
| population_selection_association | provenance_non_penalty |
| possible_genetic_background_confound | caveat_penalty |
| pre_infection_baseline | provenance_non_penalty |
| qc_reviewed | provenance_non_penalty |
| quantseq_3prime_tag_seq | provenance_non_penalty |
| raw_reanalysis | provenance_non_penalty |
| replicate_culture_independence_unverified | caveat_penalty |
| selected_A_pc2_displacement | caveat_penalty |
| sensitivity_direction_robust | provenance_non_penalty |
| sequencing_pool_partly_confounded | caveat_penalty |
| shared_control_group | caveat_penalty |
| shared_control_group_with_PRJEB86618 | caveat_penalty |
| stage_blocked_contrast | provenance_non_penalty |
| tank_random_effect_mixed_model | provenance_non_penalty |
| temperature_response | provenance_non_penalty |
| time_averaged_over_diurnal_cycle | provenance_non_penalty |
| timepoint_blocked_contrast | provenance_non_penalty |
| tire_leachate | provenance_non_penalty |
| tire_microparticles | provenance_non_penalty |
| unpaired_primary_tank_paired_sensitivity | provenance_non_penalty |
| vibrio_12h | provenance_non_penalty |
| vibrio_24h | provenance_non_penalty |
| vibrio_species_not_stated | caveat_penalty |
| virulent_vibrio_crassostreae | provenance_non_penalty |
| virulent_vibrio_tasmaniensis | provenance_non_penalty |
