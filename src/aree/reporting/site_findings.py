"""Study-level website findings derived only from committed real evidence.

No effects are pooled across studies and no cross-context candidate score is used.
"""
import pandas as pd

from aree.io import read_tsv
from aree.meta_analysis.random_effects import inference_eligible
from aree.reporting.evidence_cards import card_filename

FINDINGS_PER_STUDY = 10
EVIDENCE_COLUMNS = [
    'study_id', 'feature_id_standardized', 'effect_size', 'effect_size_type',
    'standard_error', 'p_value', 'adjusted_p_value', 'sample_size', 'sample_comparison',
    'phenotype', 'stressor', 'tissue', 'life_stage', 'resilience_classification',
    'mapping_confidence', 'quality_flags', 'source_file', 'input_checksum',
    'genome_assembly', 'annotation_version', 'workflow_version',
]


def load_findings(root, studies):
    """Return all real evidence, study summaries and a reproducible featured subset."""
    tables, summaries, featured = [], [], []
    for study_id, study in sorted(studies.items()):
        if study['data_availability']['status'].startswith('simulated'):
            continue
        path = root / 'data/harmonized' / (study_id + '.tsv')
        if not path.exists():
            summaries.append(dict(study_id=study_id, records=0, tested=0, significant=None,
                                  higher=None, lower=None, status='Not harmonized; see design limitations'))
            continue
        table = read_tsv(path, usecols=EVIDENCE_COLUMNS)
        if not table.study_id.eq(study_id).all() or not table.sample_size.eq(study['sample_size']).all():
            raise ValueError('Study ID or final sample count does not match registry: ' + study_id)
        if table.feature_id_standardized.duplicated().any():
            raise ValueError('Multiple contrasts per feature require explicit selection: ' + study_id)
        annotations = root / 'reports/analysis' / study['accessions']['bioproject'] / 'refseq_gene_annotations.tsv'
        table['description'] = 'Annotation unavailable'
        if annotations.exists():
            annotation = read_tsv(annotations, usecols=['feature_id_standardized', 'description'])
            if annotation.feature_id_standardized.duplicated().any():
                raise ValueError('Ambiguous gene annotations: ' + study_id)
            descriptions = annotation.set_index('feature_id_standardized').description
            table['description'] = table.feature_id_standardized.map(descriptions).fillna('Annotation unavailable')
        eligible = inference_eligible(table)
        table['inference_status'] = eligible.map({True: 'Study-level association; see design caveats', False: 'Exploratory; tank clustering unmodeled'})
        significant = table[eligible & table.adjusted_p_value.lt(0.05)]
        exploratory = not eligible.any()
        summaries.append(dict(study_id=study_id, records=len(table), tested=int(table.p_value.notna().sum()),
                              significant=None if exploratory else len(significant),
                              higher=None if exploratory else int(significant.effect_size.gt(0).sum()),
                              lower=None if exploratory else int(significant.effect_size.lt(0).sum()),
                              status='Exploratory; nominal FDR not used' if exploratory else 'Study-level associations'))
        table['selection_basis'] = 'FDR < 0.05; ordered by adjusted p-value, then absolute effect'
        if exploratory:
            selected = table.assign(_magnitude=table.effect_size.abs()).sort_values(
                ['_magnitude', 'feature_id_standardized'], ascending=[False, True]).head(FINDINGS_PER_STUDY).copy()
            selected['selection_basis'] = 'Descriptive absolute effect; nominal FDR not used for selection'
        else:
            selected = significant.assign(_magnitude=significant.effect_size.abs()).sort_values(
                ['adjusted_p_value', '_magnitude', 'feature_id_standardized'], ascending=[True, False, True]
            ).head(FINDINGS_PER_STUDY).copy()
        tables.append(table)
        featured.append(selected.drop(columns=['_magnitude'], errors='ignore'))
    known_files = {p.stem for p in (root / 'data/harmonized').glob('*.tsv')}
    real_ids = {s for s in studies if not studies[s]['data_availability']['status'].startswith('simulated')}
    if known_files - real_ids:
        raise ValueError('Unknown or simulated study in real evidence directory')
    evidence = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=EVIDENCE_COLUMNS + ['description', 'inference_status', 'selection_basis'])
    highlights = pd.concat(featured, ignore_index=True) if featured else evidence.copy()
    return evidence, pd.DataFrame(summaries), highlights


def number(value):
    return 'not available' if pd.isna(value) else '{:.4g}'.format(float(value))


def findings_table(highlights, md_table, prefix=''):
    rows = []
    for row in highlights.itertuples(index=False):
        rows.append([
            '[{}]({}evidence_cards/{})'.format(row.feature_id_standardized, prefix, card_filename(row.feature_id_standardized)),
            row.description, row.effect_size_type.replace('_', ' '), number(row.effect_size), number(row.standard_error),
            number(row.adjusted_p_value) + (' (nominal only)' if row.inference_status.startswith('Exploratory') else ''),
            row.inference_status,
        ])
    return md_table(pd.DataFrame(rows, columns=['Feature / card', 'RefSeq description', 'Effect scale', 'Effect', 'Standard error', 'Adjusted p-value', 'Interpretation']))


def real_card(feature, rows, studies, md_table):
    parts = ['---\ntitle: "Evidence card: {}"\n---\n'.format(feature),
             '**Real public-study evidence. Association evidence only; not a validated biomarker.**\n',
             '[How to interpret expression as evidence for resilience](../docs/interpreting-evidence.md#how-expression-findings-relate-to-resilience). Differential expression alone does not establish a resilience marker; higher expression does not necessarily mean greater resilience.\n',
             'This feature was selected among ten displayed findings in at least one study. All available rows for it are shown below, including nonsignificant and exploratory results. '
             'Presence in multiple datasets does not establish independent replication: exposures, tissues, life stages, and phenotype meanings differ. No cross-study score or pooled effect is calculated.\n',
             '## Study-specific effects\n']
    effect_rows = []
    for row in rows.itertuples(index=False):
        effect_rows.append(['[{}](../studies/{}.qmd)'.format(studies[row.study_id]['accessions']['bioproject'], row.study_id), row.description,
                            row.sample_comparison, row.effect_size_type.replace('_', ' '), number(row.effect_size),
                            number(row.standard_error),
                            number(row.adjusted_p_value) + (' (nominal only)' if row.inference_status.startswith('Exploratory') else ''),
                            row.inference_status])
    parts.append(md_table(pd.DataFrame(effect_rows, columns=['Study / design', 'RefSeq description', 'Contrast', 'Effect scale', 'Effect', 'Standard error', 'Adjusted p-value', 'Interpretation'])))
    for row in rows.itertuples(index=False):
        study = studies[row.study_id]
        parts.extend(['\n## {}\n'.format(row.study_id),
                      '**Phenotype interpretation:** ' + study['phenotype_direction'] + '\n',
                      '**Biological replication:** ' + study['biological_replication'] + '\n',
                      '**Evidence classification:** {}. Tissue: {}; life stage: {}.\n'.format(row.resilience_classification, row.tissue, row.life_stage),
                      '**Quality flags:** `{}`.\n'.format(row.quality_flags)])
        if row.inference_status.startswith('Exploratory'):
            parts.append('**Exploratory uncertainty:** Tank clustering is unmodeled. Standard errors and adjusted p-values above are nominal only; this row is excluded from pooling and significance rewards.\n')
        parts.append('**Study limitations:**\n')
        parts.extend('- ' + text for text in study['limitations'])
        parts.extend(['\n**Provenance:** Reference `{}` / `{}`; mapping `{}`; workflow `{}`.\n'.format(row.genome_assembly, row.annotation_version, row.mapping_confidence, row.workflow_version),
                      '[Harmonized evidence](https://github.com/sr320/AREE-co/raw/main/data/harmonized/{}.tsv). '.format(row.study_id) +
                      '[Processed input](https://github.com/sr320/AREE-co/blob/main/{}). Input SHA-256: `{}`.\n'.format(row.source_file, row.input_checksum)])
    parts.extend(['\n## Next validation step\n', 'Independent biological replication with matched exposure and measured phenotype, an appropriate experimental-unit model, and targeted validation in the relevant tissue and life stage.\n'])
    return '\n'.join(parts)


def write_real_findings(out, evidence, summaries, highlights, studies, md_table, write):
    """Publish a bounded browsing set plus full study-level downloads."""
    features = sorted(set(highlights.feature_id_standardized))
    for feature, rows in evidence[evidence.feature_id_standardized.isin(features)].groupby('feature_id_standardized'):
        write(out / 'evidence_cards' / card_filename(feature), real_card(feature, rows, studies, md_table))
    cards = pd.DataFrame([['[{}]({})'.format(feature, card_filename(feature)),
                           ', '.join(sorted(set(highlights.loc[highlights.feature_id_standardized == feature, 'study_id'])))]
                          for feature in features], columns=['Real evidence card', 'Selected in studies'])
    write(out / 'evidence_cards/index.qmd', '---\ntitle: "Real-study evidence cards"\n---\n\n'
          'Cards cover the union of ten displayed findings per harmonized study; each retains all available study rows and design limitations. '
          'They are selected associations, not a complete gene catalog or validated biomarkers.\n\n' + md_table(cards) + '\n')
    parts = ['---\ntitle: "Real-study findings"\n---\n',
             'Findings are reported within each study. For inference-eligible studies, ten genes are selected from adjusted p < 0.05, ordered by adjusted p-value, then absolute effect and feature ID. '
             'For exploratory PRJNA735889, ten descriptive effects are ordered by absolute effect; nominal FDR is not used for selection. '
             'These are browsing selections, not global biomarker rankings. Full harmonized tables are linked for every available study. '
             'No cross-study effect pooling is presented because matching effect scales and broad phenotype labels alone do not establish comparable contrasts.\n',
             '[How do expression findings relate to resilience?](docs/interpreting-evidence.md#how-expression-findings-relate-to-resilience) Differential expression identifies molecular differences; linking them to survival, performance or recovery requires outcome evidence and independent validation. Higher expression does not necessarily mean greater resilience.\n',
             '[Download displayed findings](downloads/real_findings.tsv) · [Download study summary](downloads/real_study_summary.tsv)\n']
    for row in summaries.itertuples(index=False):
        study = studies[row.study_id]
        parts.extend(['## {}\n'.format(study['accessions']['bioproject']),
                      '[Study design, replication and limitations](studies/{}.qmd)\n'.format(row.study_id),
                      study['phenotype_direction'] + '\n'])
        if row.records == 0:
            parts.append('**Not harmonized.** Study-level descriptive work exists, but the design does not support the released inferential evidence table. See the study limitations.\n')
            continue
        parts.append('Records: {:,}; genes with nonmissing raw p-values: {:,}. {}.\n'.format(row.records, row.tested, row.status))
        if pd.notna(row.significant):
            parts.append('Adjusted p < 0.05: {:,} ({:,} positive, {:,} negative). Counts are study-level molecular associations.\n'.format(int(row.significant), int(row.higher), int(row.lower)))
        else:
            parts.append('**Exploratory:** Tank clustering is unmodeled. Nominal significance and uncertainty may overstate precision; descriptive findings are shown without inferential claims.\n')
        parts.append('**Biological replication:** ' + study['biological_replication'] + '\n')
        parts.append('**Treatment:** ' + study['experimental_treatment'] + '\n')
        parts.append('**Control:** ' + study['control_condition'] + '\n')
        parts.append('Positive effects indicate higher expression in the treatment relative to the control; negative effects indicate lower expression.\n')
        parts.append('[Download all harmonized evidence](https://github.com/sr320/AREE-co/raw/main/data/harmonized/{}.tsv)\n'.format(row.study_id))
        parts.append(findings_table(highlights[highlights.study_id == row.study_id], md_table) + '\n')
    write(out / 'candidates.qmd', '\n'.join(parts))
    (out / 'downloads').mkdir(exist_ok=True)
    highlights.to_csv(out / 'downloads/real_findings.tsv', sep='\t', index=False)
    summaries.to_csv(out / 'downloads/real_study_summary.tsv', sep='\t', index=False)
    return features
