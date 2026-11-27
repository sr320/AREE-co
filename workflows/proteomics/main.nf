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

process PROTEOMICS_LIMMA {
    label 'bio'
    publishDir "${params.outdir}/differential", mode: 'copy'

    input:
    path abundance_table
    path design

    output:
    path "proteomics/*_differential_proteins.tsv", emit: results
    path "proteomics/analysis_method.txt", emit: analysis_method

    script:
    def log2Arg = params.log2_transformed ? "--log2-transformed=yes" : "--log2-transformed=no"
    def covArg = params.covariate ? "--covariate=${params.covariate}" : ""
    """
    mkdir -p proteomics
    Rscript ${projectDir}/scripts/run_proteomics_limma.R \\
        ${abundance_table} ${design} proteomics \\
        --reference=${params.reference_level} \\
        --test=${params.test_level} \\
        --min-samples=${params.min_samples} \\
        ${log2Arg} \\
        ${covArg}
    """
}

process EXPORT_PROTEOMICS_EVIDENCE {
    label 'aree'
    publishDir "${params.outdir}/evidence", mode: 'copy'

    input:
    path results
    path analysis_method

    output:
    path "${params.study_id}_proteomics.tsv", emit: processed
    path "${params.study_id}_mapping.tsv", emit: mapping
    path "${params.study_id}_annotations.tsv", emit: annotations

    script:
    def flags = params.quality_flags.tokenize(',').collect { "--quality-flag ${shQuote(it.trim())}" }.join(' ')
    def comparison = params.sample_comparison ?: "${params.test_level}_vs_${params.reference_level}"
    """
    aree export-assay-evidence \\
        --results ${results} \\
        --processed ${params.study_id}_proteomics.tsv \\
        --mapping ${params.study_id}_mapping.tsv \\
        --annotations ${params.study_id}_annotations.tsv \\
        --sample-comparison ${shQuote(comparison)} \\
        --feature-type protein \\
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

    PROTEOMICS_LIMMA(abund, design)
    EXPORT_PROTEOMICS_EVIDENCE(PROTEOMICS_LIMMA.out.results, PROTEOMICS_LIMMA.out.analysis_method)
    if (params.validate_evidence) {
        VALIDATE_EVIDENCE(
            EXPORT_PROTEOMICS_EVIDENCE.out.processed,
            EXPORT_PROTEOMICS_EVIDENCE.out.mapping,
            EXPORT_PROTEOMICS_EVIDENCE.out.annotations,
        )
    }
}
