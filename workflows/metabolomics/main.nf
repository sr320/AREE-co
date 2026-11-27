nextflow.enable.dsl = 2

include { VALIDATE_EVIDENCE } from '../../modules/validate_evidence'

def required(name) {
    if (!params[name]) {
        error "Missing required parameter --${name}"
    }
    return params[name]
}

def shQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process METABOLOMICS_DIFF {
    label 'bio'
    publishDir "${params.outdir}/differential", mode: 'copy'

    input:
    path abundance_table
    path design

    output:
    path "metabolomics/*_differential_metabolites.tsv", emit: results
    path "metabolomics/analysis_method.txt", emit: analysis_method

    script:
    def covArg = params.covariate ? "--covariate=${params.covariate}" : ""
    """
    mkdir -p metabolomics
    Rscript ${projectDir}/scripts/run_metabolomics_diff.R \\
        ${abundance_table} ${design} metabolomics \\
        --reference=${params.reference_level} \\
        --test=${params.test_level} \\
        --min-samples=${params.min_samples} \\
        ${covArg}
    """
}

process EXPORT_METABOLOMICS_EVIDENCE {
    label 'aree'
    publishDir "${params.outdir}/evidence", mode: 'copy'

    input:
    path results
    path analysis_method

    output:
    path "${params.study_id}_metabolomics.tsv", emit: processed
    path "${params.study_id}_mapping.tsv", emit: mapping
    path "${params.study_id}_annotations.tsv", emit: annotations

    script:
    def flags = params.quality_flags.tokenize(',').collect { "--quality-flag ${shQuote(it.trim())}" }.join(' ')
    def comparison = params.sample_comparison ?: "${params.test_level}_vs_${params.reference_level}"
    """
    aree export-assay-evidence \\
        --results ${results} \\
        --processed ${params.study_id}_metabolomics.tsv \\
        --mapping ${params.study_id}_mapping.tsv \\
        --annotations ${params.study_id}_annotations.tsv \\
        --sample-comparison ${shQuote(comparison)} \\
        --feature-type metabolite \\
        --effect-size-type log2_fold_change \\
        --direction-mode abundance \\
        --reference-release ${shQuote(params.reference_release)} \\
        --analysis-method "\$(cat ${analysis_method})" \\
        ${flags}
    """
}

workflow {
    required('study_id')
    required('test_level')
    abund = file(required('abundance_table'), checkIfExists: true)
    design = file(required('design'), checkIfExists: true)

    METABOLOMICS_DIFF(abund, design)
    EXPORT_METABOLOMICS_EVIDENCE(METABOLOMICS_DIFF.out.results, METABOLOMICS_DIFF.out.analysis_method)
    if (params.validate_evidence) {
        VALIDATE_EVIDENCE(
            EXPORT_METABOLOMICS_EVIDENCE.out.processed,
            EXPORT_METABOLOMICS_EVIDENCE.out.mapping,
            EXPORT_METABOLOMICS_EVIDENCE.out.annotations,
        )
    }
}
