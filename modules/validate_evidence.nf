process VALIDATE_EVIDENCE {
    tag "${params.study_id}"
    label 'aree'
    publishDir "${params.outdir}/evidence", mode: 'copy'

    input:
    path processed
    path mapping
    path annotation

    output:
    path "validation_report.json", emit: report

    script:
    """
    aree validate-processed \\
        --processed ${processed} \\
        --mapping ${mapping} \\
        --annotation ${annotation} \\
        --out validation_report.json
    """

    stub:
    """
    echo '{"status": "stub"}' > validation_report.json
    """
}
