from pathlib import Path

import pandas as pd
import pytest
import yaml

from aree.reporting.site_findings import load_findings, real_card
from scripts.build_site import md_table, load_registry

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def findings():
    studies = {s['study_id']: s for s in (yaml.safe_load(p.read_text()) for p in (ROOT / 'registry/studies').glob('*.yaml'))}
    return studies, load_findings(ROOT, studies)


def test_real_findings_do_not_use_simulated_data_or_global_rankings(findings):
    studies, (evidence, summaries, selected) = findings
    assert len(evidence) == 364777
    assert len(summaries) == 18
    assert len(selected) == 168
    assert selected.study_id.nunique() == 17
    # Each study features ten findings, or all of them when fewer reach FDR < 0.05.
    sizes = selected.groupby('study_id').size()
    assert sizes.le(10).all()
    assert sizes.drop('CGIG_TIREPARTICLE_RNASEQ_PRJNA856813').eq(10).all()
    assert sizes['CGIG_TIREPARTICLE_RNASEQ_PRJNA856813'] == 8
    assert all(not studies[study]['data_availability']['status'].startswith('simulated') for study in evidence.study_id.unique())
    assert 'score' not in selected
    assert 'CGIG_OA_RNASEQ_PRJNA1196326' not in set(evidence.study_id)


def test_significance_counts_and_selection_match_underlying_evidence(findings):
    _, (evidence, summaries, selected) = findings
    for row in summaries.itertuples(index=False):
        group = evidence[evidence.study_id == row.study_id]
        assert row.records == len(group)
        if group.empty or pd.isna(row.significant):
            continue
        significant = group[group.adjusted_p_value < 0.05]
        assert row.significant == len(significant)
        assert row.higher == significant.effect_size.gt(0).sum()
        assert row.lower == significant.effect_size.lt(0).sum()
        top = selected[selected.study_id == row.study_id]
        assert top.adjusted_p_value.lt(0.05).all()
        assert top.adjusted_p_value.is_monotonic_increasing
        assert set(top.feature_id_standardized).issubset(set(significant.feature_id_standardized))


def test_exploratory_selection_ignores_nominal_significance(findings):
    _, (evidence, summaries, selected) = findings
    study = 'CGIG_OA_RNASEQ_PRJNA735889'
    row = summaries.set_index('study_id').loc[study]
    assert pd.isna(row.significant)
    top = selected[selected.study_id == study]
    original = evidence[evidence.study_id == study].assign(magnitude=lambda x: x.effect_size.abs())
    expected = original.sort_values(['magnitude', 'feature_id_standardized'], ascending=[False, True]).head(10)
    assert list(top.feature_id_standardized) == list(expected.feature_id_standardized)
    assert top.selection_basis.str.contains('nominal FDR not used').all()


def test_real_card_preserves_all_observations_and_shared_control_caveats(findings):
    studies, (evidence, _, selected) = findings
    feature = 'NCBI:GeneID:105327837'
    assert feature in set(selected.feature_id_standardized)
    rows = evidence[evidence.feature_id_standardized == feature]
    card = real_card(feature, rows, studies, md_table)
    for study in rows.study_id:
        assert study in card
        assert studies[study]['biological_replication'] in card
    assert 'share a reference group' in card
    assert 'No cross-study score or pooled effect' in card
    assert 'SHA-256' in card
    assert 'nominal only' in card
    assert len(rows) == len(evidence[evidence.feature_id_standardized == feature])


def test_build_rejects_stale_sample_metadata(findings):
    studies, _ = findings
    studies = dict(studies)
    study = 'CGIG_HAB_RNASEQ_PRJEB18614'
    studies[study] = dict(studies[study], sample_size=26)
    with pytest.raises(ValueError, match='final sample count'):
        load_findings(ROOT, studies)


def test_site_rejects_stale_registry_csv(tmp_path, monkeypatch):
    import shutil
    import scripts.build_site as site

    shutil.copytree(ROOT / 'registry', tmp_path / 'registry')
    path = tmp_path / 'registry/study_registry.csv'
    registry = pd.read_csv(path)
    registry.loc[0, 'quality_control_status'] = 'stale metadata'
    registry.to_csv(path, index=False)
    monkeypatch.setattr(site, 'ROOT', tmp_path)
    with pytest.raises(ValueError, match='Registry CSV is stale'):
        load_registry()
