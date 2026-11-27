from pathlib import Path

import pandas as pd
import pytest
import yaml

from aree.io import read_tsv
from aree.meta_analysis.random_effects import EXPLORATORY_FLAG, meta_analysis_table
from aree.prioritize.scoring import score_table
from scripts.build_site import build

ROOT = Path(__file__).resolve().parents[1]


def evidence(rows):
    base = dict(feature_id_standardized='G1', feature_type='gene', effect_size_type='log2_fold_change',
                effect_size=1.0, standard_error=0.2, sample_size=12,
                resilience_classification='resilience_associated', mapping_confidence='exact',
                quality_flags='none', adjusted_p_value=0.001, tissue='gill', life_stage='adult',
                phenotype='thermal_tolerance', stressor='temperature')
    return pd.DataFrame([dict(base, **row) for row in rows])


def test_shared_controls_do_not_establish_independent_replication():
    flags = 'dependence_group=shared_controls'
    shared = score_table(evidence([dict(study_id='A', quality_flags=flags),
                                   dict(study_id='B', quality_flags=flags)])).iloc[0]
    assert shared.n_studies == 2
    assert shared.n_independent_studies == 1
    assert shared.total_biological_sample_size == 12
    assert 'lower bound' in shared.sample_size_status
    assert pd.isna(shared.direction_consistency)
    assert shared.category == 'Emerging candidate requiring replication'
    independent = score_table(evidence([dict(study_id='A'), dict(study_id='B')])).iloc[0]
    assert independent.n_independent_studies == 2
    assert independent.total_biological_sample_size == 24
    assert independent.category == 'High-priority cross-study candidate'
    assert shared.score < independent.score


def test_shared_group_plus_independent_replication_and_conflicting_directions():
    rows = [dict(study_id='A', quality_flags='dependence_group=shared'),
            dict(study_id='B', quality_flags='dependence_group=shared'),
            dict(study_id='C')]
    row = score_table(evidence(rows)).iloc[0]
    assert row.n_independent_studies == 2
    assert row.total_biological_sample_size == 24
    assert row.category == 'High-priority cross-study candidate'
    rows[-1]['effect_size'] = -1.0
    row = score_table(evidence(rows)).iloc[0]
    assert row.direction_consistency == 0.5
    assert row.category != 'High-priority cross-study candidate'


def test_inconsistent_dependence_metadata_fails_closed():
    with pytest.raises(ValueError, match='Inconsistent dependence'):
        score_table(evidence([dict(study_id='A'), dict(study_id='A', quality_flags='dependence_group=shared')]))


@pytest.mark.parametrize('third', [[], [dict(study_id='C')]])
def test_correlated_effects_are_visible_without_independence_based_pooling(third):
    from aree.meta_analysis.random_effects import random_effects

    frame = evidence([dict(study_id='A', quality_flags='dependence_group=shared'),
                      dict(study_id='B', quality_flags='dependence_group=shared')] + third)
    for result in [random_effects(frame), meta_analysis_table(frame).iloc[0]]:
        assert result['pooling_status'] == 'dependent_effects_require_covariance'
        assert pd.isna(result['pooled_effect']) and pd.isna(result['i2_percent'])
        assert result['n_effects_excluded'] == len(frame)
        assert set(result['excluded_study_ids'].split(';')) == set(frame.study_id)
    assert not score_table(frame).empty


def test_excluded_exploratory_effect_does_not_block_independent_pooling():
    frame = evidence([dict(study_id='A', quality_flags='dependence_group=shared;' + EXPLORATORY_FLAG),
                      dict(study_id='B', quality_flags='dependence_group=shared'), dict(study_id='C')])
    result = meta_analysis_table(frame).iloc[0]
    assert result.pooling_status == 'pooled'
    assert result.n_effects == 2 and result.excluded_study_ids == 'A'


def test_inconsistent_study_sample_counts_fail_before_scoring():
    with pytest.raises(ValueError, match='Inconsistent sample sizes'):
        score_table(evidence([dict(study_id='A', sample_size=24),
                              dict(study_id='A', sample_size=26, feature_id_standardized='G2')]))


def test_all_committed_evidence_counts_match_final_registry():
    registry = pd.read_csv(ROOT / 'registry/study_registry.csv').set_index('study_id')
    paths = list((ROOT / 'data/harmonized').glob('*.tsv')) + [ROOT / 'data/demo/harmonized_evidence.tsv']
    for path in paths:
        frame = read_tsv(path, usecols=['study_id', 'sample_size'])
        for study_id, group in frame.groupby('study_id'):
            study = yaml.safe_load((ROOT / 'registry/studies' / (study_id + '.yaml')).read_text())
            expected = study['sample_size']
            assert registry.loc[study_id, 'sample_size'] == expected, study_id
            assert group.sample_size.eq(expected).all(), '{} has stale sample counts in {}'.format(study_id, path)


def test_primary_counts_match_actual_manifest_exclusions():
    from scripts.prepare_prjna735889_fastqs import EXCLUDED_RUNS

    ph = read_tsv(ROOT / 'data/manifests/CGIG_OA_RNASEQ_PRJNA735889_runs.tsv')
    ph = ph[~ph.run_accession.isin(EXCLUDED_RUNS)]
    hab = read_tsv(ROOT / 'data/manifests/CGIG_HAB_RNASEQ_PRJEB18614_runs.tsv')
    hab = hab[hab.timepoint != 'd2_zt08']
    for study_id, included, per_arm in [('CGIG_OA_RNASEQ_PRJNA735889', ph, 25),
                                      ('CGIG_HAB_RNASEQ_PRJEB18614', hab, 12)]:
        study = yaml.safe_load((ROOT / 'registry/studies' / (study_id + '.yaml')).read_text())
        assert len(included) == study['sample_size'] == 2 * per_arm
        assert included.condition.value_counts().eq(per_arm).all()


def test_exploratory_significance_is_not_rewarded_or_pooled():
    rows = [dict(study_id='A', quality_flags=EXPLORATORY_FLAG, adjusted_p_value=1e-100),
            dict(study_id='B', quality_flags=EXPLORATORY_FLAG, adjusted_p_value=1e-100)]
    frame = evidence(rows)
    score = score_table(frame).iloc[0]
    assert pd.isna(score.best_adjusted_p_value)
    assert score.category != 'High-priority cross-study candidate'
    frame['adjusted_p_value'] = 0.9
    assert score_table(frame).iloc[0].score == score.score
    pooled = meta_analysis_table(frame).iloc[0]
    assert pooled.n_effects == 0
    assert pooled.n_effects_excluded == 2
    assert pooled.pooling_status == 'no_inference_eligible_effects'
    assert pd.isna(pooled.pooled_effect)
    frame.loc[1, 'quality_flags'] = 'none'
    pooled = meta_analysis_table(frame).iloc[0]
    assert pooled.n_effects == 1
    assert pooled.study_ids == 'B'


def test_committed_shared_controls_and_tank_flags_are_enforced():
    heat = read_tsv(ROOT / 'data/harmonized/CGIG_HEAT_RNASEQ_PRJEB86646.tsv')
    starvation = read_tsv(ROOT / 'data/harmonized/CGIG_STARVATION_RNASEQ_PRJEB86618.tsv')
    feature = heat.iloc[0].feature_id_standardized
    matched = pd.concat([heat[heat.feature_id_standardized == feature],
                         starvation[starvation.feature_id_standardized == feature]])
    score = score_table(matched).iloc[0]
    assert score.n_studies == 2 and score.n_independent_studies == 1
    assert score.category != 'High-priority cross-study candidate'
    ph = read_tsv(ROOT / 'data/harmonized/CGIG_OA_RNASEQ_PRJNA735889.tsv')
    assert ph.sample_size.eq(50).all()
    assert ph.quality_flags.str.contains('tank_random_effect_mixed_model').all()
    assert not ph.quality_flags.str.contains(EXPLORATORY_FLAG).any()
    assert pd.notna(score_table(ph.head(1)).iloc[0].best_adjusted_p_value)


def test_website_distinguishes_synthetic_results_and_exposes_limitations(tmp_path):
    out = build(tmp_path / 'site')
    index = (out / 'index.qmd').read_text()
    assert 'title: "Aquaculture Resilience Evidence Companion"' in index
    real_count = sum(len(read_tsv(path, usecols=['study_id'])) for path in (ROOT / 'data/harmonized').glob('*.tsv'))
    assert '>{:,}</div><div class="stat-label">Real harmonized evidence records'.format(real_count) in index
    assert '## What the work involves' in index
    assert 'Records are neither unique genes nor validated biomarkers' in index
    assert '## Top candidates' not in index
    assert '## One displayed finding per harmonized study' in index
    assert '## Simulated demo candidates' not in index
    assert '(demo/index.qmd)' in index
    studies = [yaml.safe_load(path.read_text()) for path in (ROOT / 'registry/studies').glob('*.yaml')]
    real_studies = [study for study in studies if not study['data_availability']['status'].startswith('simulated')]
    projects = {study['accessions']['bioproject'] for study in real_studies}
    assert '>{}</div><div class="stat-label">Real public BioProjects'.format(len(projects)) in index
    assert '>{}</div><div class="stat-label">Registered real-study contrasts'.format(len(real_studies)) in index
    findings = (out / 'candidates.qmd').read_text()
    assert '{#cgig_heat_rnaseq_prjna913164}' in findings
    assert '{#cgig_multistress_rnaseq_prjna913164}' in findings
    for path in [out / 'demo/candidates.qmd', out / 'demo/meta-analysis.qmd', out / 'demo/evidence_cards/index.qmd']:
        assert 'simulated' in path.read_text().lower()
    for path in (out / 'demo/evidence_cards').glob('*.md'):
        assert 'Simulated demo evidence' in path.read_text()
    assert (out / 'evidence_cards/index.qmd').exists()
    assert 'up to ten' in (out / 'evidence_cards/index.qmd').read_text()
    assert 'up to ten' in (out / 'candidates.qmd').read_text()
    assert '## Reading the results' in (out / 'docs/methods.md').read_text()
    for path in (out / 'evidence_cards').glob('*.md'):
        assert 'Real public-study evidence' in path.read_text()
        assert 'Simulated demo evidence' not in path.read_text()
    assert 'CGIG_HEAT_RNASEQ_001' not in (out / 'studies.qmd').read_text()
    ph = (out / 'studies/CGIG_OA_RNASEQ_PRJNA735889.qmd').read_text()
    assert 'tank clustering is not modeled' not in ph
    assert 'variancePartition::dream' in ph
    synth = (out / 'meta-analysis.qmd').read_text()
    assert 'acute_heat_response_gill' in synth
    assert 'not validated biomarkers' in synth
    heat = (out / 'studies/CGIG_HEAT_RNASEQ_PRJEB86646.qmd').read_text()
    assert 'share a reference group' in heat
    hab = (out / 'studies/CGIG_HAB_RNASEQ_PRJEB18614.qmd').read_text()
    assert '**Primary-analysis sample count:** 24.' in hab
    assert 'not individual oysters or independent tanks' in hab
    # The query page ships only alongside its data, and keeps the interpretive caveats.
    config = (out / '_quarto.yml').read_text()
    assert 'title: AREC' in config
    assert 'Aquaculture Resilience Evidence Companion' in config
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        assert 'query.qmd' not in config
        assert not (out / 'query.qmd').exists()
    else:
        assert 'query.qmd' in config
        assert (out / 'downloads/aree_evidence.parquet').exists()
        assert (out / 'query.js').exists()
        query = (out / 'query.qmd').read_text()
        assert 'not validated biomarkers' in query
        assert 'nominal' in query
        assert 'id="aree-inference"' in query and 'Inference-eligible only' in query
        assert '<th>Scale</th><th>Effect</th>' in query
        assert 'Effect (log2FC)' not in query
        query_js = (out / 'query.js').read_text()
        assert 'effect_size_type' in query_js
        assert 'contains(lower(' in query_js and 'ILIKE' not in query_js
        assert "inference_status NOT LIKE 'Exploratory%'" in query_js
        assert "CASE WHEN inference_status LIKE 'Exploratory%' THEN 1 ELSE 0 END" in query_js


def test_exploratory_study_does_not_veto_qualified_candidate_or_inflate_replication():
    qualified = [
        dict(study_id='A', sample_size=20, effect_size=1.2, adjusted_p_value=1e-4),
        dict(study_id='B', sample_size=20, effect_size=1.4, adjusted_p_value=1e-4),
    ]
    exp = dict(
        study_id='EXP',
        sample_size=50,
        effect_size=1.0,
        adjusted_p_value=0.9,
        resilience_classification='environmental_exposure',
        quality_flags=EXPLORATORY_FLAG,
    )
    row = score_table(evidence(qualified + [exp])).iloc[0]
    assert row.n_studies == 3 and row.n_independent_studies == 3
    assert row.category == 'High-priority cross-study candidate'
    # Changing exploratory sample size or p-value gives no replication or significance reward.
    exp_inflated = dict(exp, sample_size=500, adjusted_p_value=1e-100)
    assert score_table(evidence(qualified + [exp_inflated])).iloc[0].score == row.score

    # One eligible study + one exploratory study cannot promote to high priority.
    one_plus_exp = score_table(evidence([qualified[0], exp])).iloc[0]
    assert pd.isna(one_plus_exp.direction_consistency)
    assert one_plus_exp.category != 'High-priority cross-study candidate'

    # Exploratory agreement cannot mask conflicting eligible studies.
    conflicting = [
        dict(study_id='A', effect_size=2.0),
        dict(study_id='B', effect_size=-2.0),
        exp,
    ]
    conflict_row = score_table(evidence(conflicting)).iloc[0]
    assert conflict_row.direction_consistency == 0.5
    assert conflict_row.category != 'High-priority cross-study candidate'


def test_provenance_quality_flags_do_not_penalize_data_quality():
    clean_none = score_table(evidence([dict(study_id='A', quality_flags='none')])).iloc[0]
    clean_provenance = score_table(
        evidence([dict(
            study_id='A',
            quality_flags='raw_reanalysis;lineage_balanced_contrast;pathogen_challenge_response;fdr_lt_0.05;dependence_group=G1',
        )])
    ).iloc[0]
    assert clean_provenance.score == clean_none.score

    for caveat in ('raw_reanalysis;pooled_libraries', 'raw_reanalysis;no_linked_publication'):
        caveated = score_table(evidence([dict(study_id='A', quality_flags=caveat)])).iloc[0]
        assert caveated.score < clean_none.score


def test_study_detail_page_tank_callout_does_not_hardcode_prjna735889_numbers():
    from scripts.build_site import study_detail_page

    study = {
        'study_id': 'CGIG_SYNTHETIC_TANK_TEST',
        'data_availability': {'status': 'public_raw_aree_reanalysis_complete'},
        'sample_size': 24,
        'analysis_status': 'raw_reanalysis_harmonized_exploratory_tank_clustering_unmodeled',
        'experimental_treatment': 'High temp',
        'control_condition': 'Control temp',
        'biological_replication': '24 oysters across 4 tanks (6 per tank)',
        'phenotype_direction': 'Response only',
        'quality_control_status': 'Passed QC',
        'limitations': ['Tank clustering unmodeled'],
        'provenance': {'source_links': ['https://example.org/study']},
    }
    page = study_detail_page(study)
    assert 'tank clustering is not modeled' in page
    assert '24 oysters across 4 tanks' in page
    assert '50 oysters' not in page and 'ten tanks' not in page
