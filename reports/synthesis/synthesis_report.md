# Real-Study Cross-Study Synthesis Report

**Association evidence only; not validated biomarkers.** Cross-study pooling is restricted to comparability-reviewed biological strata with at least two independent BioProjects and two independent replication groups, using `REML` random-effects meta-analysis and shared-control block covariance (`rho = 0.5`).

## Stratum Summary

| stratum_id | n_contrasts | n_bioprojects | n_independent_groups | multi_group_pooled_features | pooled_fdr_lt_0_05 | concordant_fdr_lt_0_05 | median_i2_percent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| acute_heat_response_gill | 4 | 3 | 3 | 22313 | 415 | 371 | 10.07 |
| thermal_protection_poms_disease_resistance | 2 | 2 | 2 | 14906 | 3899 | 3325 | 0.0 |
| vibrio_challenge_immune_response | 4 | 2 | 2 | 20242 | 818 | 816 | 0.0 |
| oshv1_challenge_immune_response | 2 | 2 | 2 | 14457 | 2705 | 2394 | 0.0 |
| ocean_acidification_whole_organism | 2 | 2 | 2 | 19179 | 16 | 16 | 0.0 |

## Estimator and Covariance Sensitivity

| stratum_id | method | shared_control_rho | n_fdr_lt_0_05 | n_concordant_fdr_lt_0_05 | median_i2_percent | pearson_r_vs_primary | jaccard_fdr_0_05_vs_primary |
| --- | --- | --- | --- | --- | --- | --- | --- |
| acute_heat_response_gill | reml | 0.5 | 415 | 371 | 10.07 | 1.0 | 1.0 |
| acute_heat_response_gill | dl | 0.5 | 382 | 361 | 12.98 | 0.9993 | 0.8155 |
| acute_heat_response_gill | reml | 0.0 | 857 | 724 | 15.12 | 0.99673 | 0.4671 |
| acute_heat_response_gill | reml | 0.8 | 341 | 306 | 9.85 | 0.98671 | 0.5429 |
| thermal_protection_poms_disease_resistance | reml | 0.5 | 3899 | 3325 | 0.0 | 1.0 | 1.0 |
| thermal_protection_poms_disease_resistance | dl | 0.5 | 3899 | 3325 | 0.0 | 1.0 | 1.0 |
| thermal_protection_poms_disease_resistance | reml | 0.0 | 3899 | 3325 | 0.0 | 1.0 | 1.0 |
| thermal_protection_poms_disease_resistance | reml | 0.8 | 3899 | 3325 | 0.0 | 1.0 | 1.0 |
| vibrio_challenge_immune_response | reml | 0.5 | 818 | 816 | 0.0 | 1.0 | 1.0 |
| vibrio_challenge_immune_response | dl | 0.5 | 818 | 816 | 0.0 | 1.0 | 1.0 |
| vibrio_challenge_immune_response | reml | 0.0 | 1282 | 1281 | 11.77 | 0.99623 | 0.6216 |
| vibrio_challenge_immune_response | reml | 0.8 | 689 | 687 | 0.0 | 0.98534 | 0.7382 |
| oshv1_challenge_immune_response | reml | 0.5 | 2705 | 2394 | 0.0 | 1.0 | 1.0 |
| oshv1_challenge_immune_response | dl | 0.5 | 2705 | 2394 | 0.0 | 1.0 | 1.0 |
| oshv1_challenge_immune_response | reml | 0.0 | 2705 | 2394 | 0.0 | 1.0 | 1.0 |
| oshv1_challenge_immune_response | reml | 0.8 | 2705 | 2394 | 0.0 | 1.0 | 1.0 |
| ocean_acidification_whole_organism | reml | 0.5 | 16 | 16 | 0.0 | 1.0 | 1.0 |
| ocean_acidification_whole_organism | dl | 0.5 | 16 | 16 | 0.0 | 1.0 | 1.0 |
| ocean_acidification_whole_organism | reml | 0.0 | 16 | 16 | 0.0 | 1.0 | 1.0 |
| ocean_acidification_whole_organism | reml | 0.8 | 16 | 16 | 0.0 | 1.0 | 1.0 |

## Acute heat response in gill (`acute_heat_response_gill`)

All four contrasts measure acute elevated-temperature exposure (30-35 C seawater vs ambient) in Crassostrea gigas gill (ctenidium) tissue across adult and juvenile oysters. PRJNA407831 6h and 24h contrasts share 0h gill pools (dependence group PRJNA407831_0h_pools) and are combined via block-covariance GLS before random-effects pooling across the three independent BioProjects.

**Limitations:**
- Excludes CGIG_THERMOTOL_RNASEQ_PRJNA694496 (constitutive larval expression after two generations of artificial selection, not acute heat shock).
- Exposure durations span 2 h to 24 h and temperatures span 30 C to 35 C; residual heterogeneity (I2) reflects exposure-time and population differences.
- Measures transcriptional heat response (stress_response), not individual survival outcomes.

### Top 10 pooled associations (concordant signals prioritized)

| feature_id_standardized | description | n_studies | n_independent_groups | pooled_effect | pooled_standard_error | adjusted_p_value | i2_percent | direction_concordant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:GeneID:105331566 | transcription factor MafK | 4 | 3 | 0.6326462569239821 | 0.04342221189376567 | 9.776561768982055e-44 | 0.0 | True |
| NCBI:GeneID:105340538 | protein hook homolog | 4 | 3 | 1.2668957182634502 | 0.11672667487355362 | 2.1407377494798957e-23 | 38.507006047252254 | True |
| NCBI:GeneID:105322994 | E3 SUMO-protein ligase EGR2 | 4 | 3 | 1.0375943457179386 | 0.10415841009848852 | 1.6675723310587705e-19 | 0.0 | True |
| NCBI:GeneID:105328771 | uncharacterized LOC105328771 | 4 | 3 | -1.0492610059586673 | 0.11749836346346071 | 2.3762701605202824e-15 | 0.0 | True |
| NCBI:GeneID:105322434 | uncharacterized LOC105322434 | 3 | 2 | 1.3799550374708385 | 0.16887835795529266 | 1.3615136246076878e-12 | 0.0 | True |
| NCBI:GeneID:105344453 | band 4.1-like protein 1 | 4 | 3 | -0.20103589575486744 | 0.024765188728012293 | 1.767103327183415e-12 | 0.0 | True |
| NCBI:GeneID:105345498 | semaphorin-1A | 4 | 3 | -0.1970885001976883 | 0.025249111160974006 | 1.8852594703458865e-11 | 0.0 | True |
| NCBI:GeneID:105333367 | uncharacterized LOC105333367 | 4 | 3 | -0.26337921322554275 | 0.034890720066988096 | 1.2262654495478557e-10 | 0.0 | True |
| NCBI:GeneID:105341538 | cyclin-dependent kinase-like 1 | 4 | 3 | -0.30284409046015726 | 0.04092691396255141 | 3.385693059935698e-10 | 0.0 | True |
| NCBI:GeneID:105322510 | innexin unc-9 | 4 | 3 | -0.19760411599738648 | 0.028319225912170162 | 6.693277556775131e-09 | 0.0 | True |

## Warm acclimation protecting against POMS / OsHV-1 (`thermal_protection_poms_disease_resistance`)

Both contrasts compare warm-acclimated oysters (29-30 C) with permissive-temperature controls (21-23 C) in whole-animal tissue from independent BioProjects where warm acclimation demonstrably reduced subsequent OsHV-1/POMS mortality (resilience_associated).

**Limitations:**
- PRJEB86646 sampled spat pools at 0 h pre-cohabitation, whereas PRJNA593309 sampled individual juveniles (n = 5) at 12 h post-cohabitation.
- With k = 2 independent studies, between-study tau2 and I2 have wide uncertainty; concordant direction is more informative than exact pooled magnitude.

### Top 10 pooled associations (concordant signals prioritized)

| feature_id_standardized | description | n_studies | n_independent_groups | pooled_effect | pooled_standard_error | adjusted_p_value | i2_percent | direction_concordant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:GeneID:105327837 | uncharacterized LOC105327837 | 2 | 2 | -3.790125045186602 | 0.18136022920679432 | 8.266390171462884e-93 | 0.0 | True |
| NCBI:GeneID:105337015 | uncharacterized LOC105337015 | 2 | 2 | 1.415145255082775 | 0.08465480278459445 | 7.387925250521933e-59 | 0.0 | True |
| NCBI:GeneID:105345603 | putative peptidyl-prolyl cis-trans isomerase | 2 | 2 | -1.8288986799066633 | 0.12160851174846335 | 2.0184346263753953e-47 | 0.0 | True |
| NCBI:GeneID:105323167 | small ribosomal subunit protein eS27-like | 2 | 2 | 0.5734421876023709 | 0.03954480988007745 | 4.444603954145585e-44 | 0.0 | True |
| NCBI:GeneID:105329390 | asparagine synthetase [glutamine-hydrolyzing] | 2 | 2 | -2.1274179047160184 | 0.14826275018806748 | 3.2235613216648405e-43 | 0.0 | True |
| NCBI:GeneID:105347892 | hillarin | 2 | 2 | 2.786515576231872 | 0.21494450295235715 | 4.870869474768184e-35 | 0.0 | True |
| NCBI:GeneID:105327885 | GTP-binding nuclear protein Ran | 2 | 2 | 0.6476660027938647 | 0.050044344775104936 | 5.563955374321165e-35 | 0.0 | True |
| NCBI:GeneID:105343175 | F-box only protein 43 | 2 | 2 | 3.2411253101850317 | 0.25731178764653667 | 4.133473909567877e-33 | 0.0 | True |
| NCBI:GeneID:105341700 | glutamine--fructose-6-phosphate aminotransferase [isomerizing] 2 | 2 | 2 | -0.872029851578285 | 0.0733278657584372 | 2.151992028891422e-29 | 0.0 | True |
| NCBI:GeneID:105325905 | developmentally-regulated GTP-binding protein 1 | 2 | 2 | 0.715791200300656 | 0.060579453626573014 | 4.824746665466134e-29 | 0.0 | True |

## Acute Vibrio challenge immune response (`vibrio_challenge_immune_response`)

Four acute Vibrio challenge contrasts across two independent BioProjects. Within each BioProject, contrasts share a control arm (PRJNA515169_mock_pools and PRJNA852380_controls) and are collapsed via shared-control block covariance before pooling across the two independent studies.

**Limitations:**
- PRJNA515169 assayed whole spat 8 h after intramuscular injection of characterized virulent strains, whereas PRJNA852380 assayed juvenile digestive gland at 12 h and 24 h with unrecorded Vibrio species/route.
- Represents pathogen-challenge response (disease_associated), not survival resistance between families.

### Top 10 pooled associations (concordant signals prioritized)

| feature_id_standardized | description | n_studies | n_independent_groups | pooled_effect | pooled_standard_error | adjusted_p_value | i2_percent | direction_concordant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:GeneID:136273704 | coadhesin-like | 4 | 2 | 3.1332106114374008 | 0.1264444574006656 | 3.030757243479142e-131 | 0.0 | True |
| NCBI:GeneID:136276560 | uncharacterized LOC136276560 | 4 | 2 | 3.8441904918212795 | 0.20987436136548598 | 6.171421116686859e-71 | 0.0 | True |
| NCBI:GeneID:105322542 | interleukin 17-like protein | 4 | 2 | 3.7115056168660656 | 0.22081094571903523 | 1.0696706055654979e-59 | 0.0 | True |
| NCBI:GeneID:105329610 | FMRFamide receptor | 4 | 2 | 2.2273306883851687 | 0.13245602157487227 | 1.0696706055654979e-59 | 0.0 | True |
| NCBI:GeneID:136276180 | uncharacterized LOC136276180 | 4 | 2 | 1.7200769902054767 | 0.11443975941905853 | 1.878839974873769e-47 | 0.0 | True |
| NCBI:GeneID:105332146 | serum amyloid A-5 protein | 4 | 2 | 2.642329855162873 | 0.18123023314168124 | 1.270955040731381e-44 | 0.0 | True |
| NCBI:GeneID:105321692 | cytochrome P450 2C31-like | 4 | 2 | 2.339426814090576 | 0.17559196483619907 | 4.914054957613966e-37 | 45.567072588916325 | True |
| NCBI:GeneID:117683663 | uncharacterized LOC117683663 | 4 | 2 | 1.0996075292420548 | 0.08314960977877175 | 1.604244677309764e-36 | 0.0 | True |
| NCBI:GeneID:105344708 | deleted in malignant brain tumors 1 protein | 4 | 2 | 1.4904853864091516 | 0.11327289685242387 | 3.426746697210382e-36 | 0.0 | True |
| NCBI:GeneID:105342093 | zinc metalloproteinase nas-13 | 4 | 2 | 1.3166013827803804 | 0.10669789039513491 | 1.1230586620470353e-31 | 0.0 | True |

## Acute OsHV-1 infection immune response (`oshv1_challenge_immune_response`)

Two independent BioProjects measuring acute whole-animal transcriptional response to OsHV-1 microvariant challenge relative to unchallenged controls in juvenile/spat Pacific oysters.

**Limitations:**
- PRJNA1329250 pools 42 spat libraries across three OsHV-1 strains and two lineages, whereas PRJNA593309 compares n = 6 juveniles (48 h vs 0 h cohabitation at 21 C).
- Reflects active viral infection response (disease_associated), not baseline genetic resistance.

### Top 10 pooled associations (concordant signals prioritized)

| feature_id_standardized | description | n_studies | n_independent_groups | pooled_effect | pooled_standard_error | adjusted_p_value | i2_percent | direction_concordant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:GeneID:105341146 | uncharacterized LOC105341146 | 2 | 2 | 3.077050020767415 | 0.23859655823587397 | 6.813251912634402e-34 | 0.0 | True |
| NCBI:GeneID:105343525 | myb-like protein X | 2 | 2 | 3.187447044609969 | 0.251336852434381 | 5.380521324792489e-33 | 0.0 | True |
| NCBI:GeneID:105336469 | NFX1-type zinc finger-containing protein 1 | 2 | 2 | 2.4824808892117667 | 0.21635689011356377 | 8.58783274302935e-27 | 0.0 | True |
| NCBI:GeneID:136269423 | uncharacterized LOC136269423 | 2 | 2 | 3.226979785642594 | 0.2905117954942618 | 3.8759750406628565e-25 | 10.739652843250234 | True |
| NCBI:GeneID:136275502 | uncharacterized LOC136275502 | 2 | 2 | 3.1445371371388213 | 0.28344314264590115 | 3.8759750406628565e-25 | 0.0 | True |
| NCBI:GeneID:105342383 | rhamnose-binding lectin | 2 | 2 | 3.1954082954190257 | 0.29130291237825806 | 1.2926285341997427e-24 | 0.0 | True |
| NCBI:GeneID:105340146 | 3'-5' exoribonuclease HELZ2 | 2 | 2 | 2.576459988555412 | 0.23704885162135994 | 3.3487775949524186e-24 | 0.0 | True |
| NCBI:GeneID:105338468 | uncharacterized LOC105338468 | 2 | 2 | 2.355737167879546 | 0.22029655020908354 | 1.9732558338373914e-23 | 0.0 | True |
| NCBI:GeneID:136276032 | uncharacterized LOC136276032 | 2 | 2 | 2.8297462426253115 | 0.26563989488305056 | 2.724753554320027e-23 | 0.0 | True |
| NCBI:GeneID:136269478 | uncharacterized LOC136269478 | 2 | 2 | 3.4111876278422453 | 0.32159763902491734 | 3.996936674554497e-23 | 0.0 | True |

## Whole-organism ocean acidification exposure (`ocean_acidification_whole_organism`)

Two independent whole-organism ocean acidification BioProjects unlocked by tank-aware mixed-model reanalysis of PRJNA735889 (~ gametogenesis + condition + (1 | tank)).

**Limitations:**
- PRJNA298285 contrasts larval pools at pH 7.9 vs 8.1 (with unlinked publication), whereas PRJNA735889 contrasts 5-month-old juveniles at pHT 6.5-6.8 vs 7.4-7.8 without an ambient pH 8.0 reference.
- Excludes CGIG_OA_RNASEQ_PRJNA826964 (adult hepatopancreas, metabolic_resilience, resequenced RNA pools).

### Top 10 pooled associations (concordant signals prioritized)

| feature_id_standardized | description | n_studies | n_independent_groups | pooled_effect | pooled_standard_error | adjusted_p_value | i2_percent | direction_concordant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NCBI:GeneID:105342907 | uncharacterized LOC105342907 | 2 | 2 | -0.7067692538857896 | 0.10936592886167144 | 1.976153146123122e-06 | 0.0 | True |
| NCBI:GeneID:105334837 | uncharacterized LOC105334837 | 2 | 2 | 0.6887465295890975 | 0.1290850736974865 | 0.0009131915159118832 | 0.0 | True |
| NCBI:GeneID:105331904 | uncharacterized LOC105331904 | 2 | 2 | 0.3659296604945886 | 0.07330782504626918 | 0.003826494014539236 | 0.0 | True |
| NCBI:GeneID:105340879 | mucin-2 | 2 | 2 | -0.5419301555922267 | 0.11244884484432414 | 0.005660312263075116 | 0.0 | True |
| NCBI:GeneID:105346061 | uncharacterized LOC105346061 | 2 | 2 | 0.49988003771172307 | 0.10382792406760785 | 0.005660312263075116 | 0.0 | True |
| NCBI:GeneID:105325527 | uncharacterized LOC105325527 | 2 | 2 | -0.4699010981143479 | 0.09838115696167969 | 0.005706433193556226 | 0.0 | True |
| NCBI:GeneID:105345814 | CCAAT/enhancer-binding protein beta | 2 | 2 | 0.33835230376105063 | 0.07131689403518121 | 0.005731149514411478 | 0.0 | True |
| NCBI:GeneID:105318255 | NF-kappa-B inhibitor epsilon-like | 2 | 2 | 0.37582753322000295 | 0.08185304421914043 | 0.009061762205452099 | 0.0 | True |
| NCBI:GeneID:105324334 | carnosine synthase 1 | 2 | 2 | 0.30970897743991427 | 0.0677290631905193 | 0.009061762205452099 | 0.0 | True |
| NCBI:GeneID:105326720 | uncharacterized LOC105326720 | 2 | 2 | 0.3065084438420973 | 0.06710319962724628 | 0.009061762205452099 | 0.0 | True |

## Cross-Context Comparison: Acute Heat Gill Meta-Analysis vs PRJNA694496 Constitutive Selection (Not Pooled)

Consistent with AREC's biological boundaries, constitutive larval expression in thermotolerance-selected oysters (`CGIG_THERMOTOL_RNASEQ_PRJNA694496`) is **not** pooled with acute adult/juvenile gill heat-shock responses (`acute_heat_response_gill`). Across 17,283 common genes: Pearson $r = 0.1034$, Spearman $\rho = 0.1313$; 182 genes reach FDR < 0.05 in both (128 same direction, 54 opposite direction).

## Unpooled Harmonized Contrasts

Ten harmonized contrasts are intentionally withheld from cross-study pooling:

- `CGIG_THERMOTOL_RNASEQ_PRJNA694496`: Constitutive ambient-temperature comparison between F2 thermotolerance-selected and control larvae; evaluated as a cross-context comparator against acute_heat_response_gill rather than pooled with acute heat shock.
- `CGIG_PESTO_F0_RNASEQ_PRJEB58194`: Same single BioProject (PRJEB58194) and laboratory pesticide mixture as CGIG_PESTO_F1DIRECT_RNASEQ_PRJEB58194; represents within-project intergenerational contrast rather than independent cross-study replication.
- `CGIG_PESTO_F1DIRECT_RNASEQ_PRJEB58194`: Same single BioProject (PRJEB58194) and laboratory pesticide mixture as CGIG_PESTO_F0_RNASEQ_PRJEB58194; lacks an independent second BioProject.
- `CGIG_TIRELEACHATE_RNASEQ_PRJNA856813`: Single BioProject and single shared-control dependence group (PRJNA856813_gill_controls) with CGIG_TIREPARTICLE_RNASEQ_PRJNA856813; cannot establish cross-study replication.
- `CGIG_TIREPARTICLE_RNASEQ_PRJNA856813`: Single BioProject and single shared-control dependence group (PRJNA856813_gill_controls) with CGIG_TIRELEACHATE_RNASEQ_PRJNA856813; cannot establish cross-study replication.
- `CGIG_STARVATION_RNASEQ_PRJEB86618`: Only one harmonized nutritional-limitation contrast.
- `CGIG_SALINITY_RNASEQ_PRJNA756710`: Only one harmonized salinity contrast.
- `CGIG_HAB_RNASEQ_PRJEB18614`: Only one harmonized harmful-algal-bloom contrast.
- `CGIG_OA_RNASEQ_PRJNA826964`: Adult hepatopancreas contrast (metabolic_resilience) with technical resequencing caveats; not biologically comparable to whole-organism OA strata.
- `CGIG_MULTISTRESS_RNASEQ_PRJNA913164`: Only one harmonized combined seawater-plus-emersion heatwave contrast.
