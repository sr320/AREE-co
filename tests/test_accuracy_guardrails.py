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
    assert ph.quality_flags.str.contains(EXPLORATORY_FLAG).all()
    assert pd.isna(score_table(ph.head(1)).iloc[0].best_adjusted_p_value)


def test_website_distinguishes_synthetic_results_and_exposes_limitations(tmp_path):
    out = build(tmp_path / 'site')
    index = (out / 'index.qmd').read_text()
    real_count = sum(len(read_tsv(path, usecols=['study_id'])) for path in (ROOT / 'data/harmonized').glob('*.tsv'))
    assert '>{:,}</div><div class="stat-label">Real harmonized evidence records'.format(real_count) in index
    assert '## Top candidates' not in index
    assert 'not findings from the real public studies' in index
    for path in [out / 'candidates.qmd', out / 'meta-analysis.qmd', out / 'evidence_cards/index.qmd']:
        assert 'simulated' in path.read_text().lower()
    for path in (out / 'evidence_cards').glob('*.md'):
        assert 'Simulated demo evidence' in path.read_text()
    ph = (out / 'studies/CGIG_OA_RNASEQ_PRJNA735889.qmd').read_text()
    assert 'tank clustering is not modeled' in ph
    assert 'ten tanks' in ph and 'no significance reward' in ph
    heat = (out / 'studies/CGIG_HEAT_RNASEQ_PRJEB86646.qmd').read_text()
    assert 'share a reference group' in heat
    hab = (out / 'studies/CGIG_HAB_RNASEQ_PRJEB18614.qmd').read_text()
    assert '**Primary-analysis sample count:** 24.' in hab
    assert 'not individual oysters or independent tanks' in hab
