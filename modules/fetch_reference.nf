process FETCH_REFERENCE {
    tag "${reference_id}"
    label 'aree'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val reference_id
    path catalog

    output:
    path "ref/transcripts.fa.gz", emit: transcripts
    path "ref/genomic.gff.gz", emit: gff
    path "ref/tx2gene.tsv", emit: tx2gene
    path "ref/reference_manifest.json", emit: manifest
    path "ref/decoys.txt", optional: true, emit: decoys

    script:
    def catFlag = catalog.name != 'NO_CATALOG' ? "--catalog ${catalog}" : ""
    def mirrorFlag = params.reference_mirror_dir ? "--mirror-dir ${params.reference_mirror_dir}" : ""
    def genomeFlag = params.include_genome_decoys ? "--include-genome" : ""
    """
    aree fetch-reference \\
        --reference "${reference_id}" \\
        --output-dir ref \\
        ${catFlag} \\
        ${mirrorFlag} \\
        ${genomeFlag}
    if [ -f ref/gentrome.fa.gz ]; then
        cp ref/gentrome.fa.gz ref/transcripts.fa.gz
    else
        cp ref/rna.fna.gz ref/transcripts.fa.gz
    fi
    """

    stub:
    """
    mkdir -p ref
    touch ref/transcripts.fa.gz ref/genomic.gff.gz ref/tx2gene.tsv ref/reference_manifest.json ref/decoys.txt
    """
}
