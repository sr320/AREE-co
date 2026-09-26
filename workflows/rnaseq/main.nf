nextflow.enable.dsl = 2

/*
 * AREE RNA-seq reanalysis: FastQC -> Salmon (decoy-aware selective alignment) -> tximport/DESeq2
 * -> AREE processed evidence, exact-ID mappings and RefSeq annotations, ready for `aree harmonize`.
 *
 * Reads are not trimmed: Salmon's selective alignment soft-clips adapters and low-quality ends,
 * matching the published AREE reanalyses of PRJNA694496 and PRJNA516762.
 */

def required(name) {
    if (!params[name]) {
        error "Missing required parameter --${name} (see workflows/rnaseq/README.md)"
    }
    return params[name]
}

process FASTQC {
    tag "${sample}"
    label 'bio'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${fastq_1} ${fastq_2}
    """

    stub:
    """
    touch ${sample}_1_fastqc.zip ${sample}_1_fastqc.html ${sample}_2_fastqc.zip ${sample}_2_fastqc.html
    """
}

process SALMON_INDEX {
    label 'bio'
    publishDir "${params.outdir}/reference", mode: 'copy', enabled: params.publish_index

    input:
    path transcripts
    path decoys

    output:
    path "salmon_index"

    script:
    def decoy_arg = decoys.name == 'NO_DECOYS' ? '' : "--decoys ${decoys}"
    """
    salmon index --transcripts ${transcripts} ${decoy_arg} --keepDuplicates \\
        --index salmon_index --threads ${task.cpus}
    """

    stub:
    """
    mkdir salmon_index && touch salmon_index/versionInfo.json
    """
}

process SALMON_QUANT {
    tag "${sample}"
    label 'bio'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)
    path index

    output:
    path "${sample}"

    script:
    """
    salmon quant --index ${index} --libType A --mates1 ${fastq_1} --mates2 ${fastq_2} \\
        --threads ${task.cpus} --validateMappings --seqBias --gcBias --output ${sample}
    """

    stub:
    """
    mkdir -p ${sample}/aux_info
    touch ${sample}/quant.sf
    echo '{"num_processed": 1, "num_mapped": 1, "percent_mapped": 100.0}' > ${sample}/aux_info/meta_info.json
    """
}

process SALMON_QC {
    label 'aree'
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'quant/*'
    path design

    output:
    path "salmon_qc_summary.tsv"

    script:
    """
    python ${params.scripts_dir}/summarize_salmon_qc.py --quant-dir quant --design-sheet ${design} \\
        --output salmon_qc_summary.tsv
    """

    stub:
    """
    touch salmon_qc_summary.tsv
    """
}

process DESEQ2 {
    label 'bio'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'quant/*'
    path design
    path tx2gene

    output:
    path "deseq2/${params.test_level}_vs_${params.reference_level}_deseq2_all_genes.tsv", emit: results
    path "deseq2/analysis_method.txt", emit: analysis_method
    path "deseq2/*", emit: all

    script:
    def ma_title = params.ma_title ? "'--ma-title=${params.ma_title}'" : ''
    """
    Rscript ${params.scripts_dir}/run_salmon_tximport_deseq2.R quant ${design} ${tx2gene} deseq2 \\
        --reference=${params.reference_level} --test=${params.test_level} \\
        --replicates=${params.replicates} --min-samples=${params.min_samples} ${ma_title}

    # Record the versions that actually ran, for the evidence rows' analysis_method.
    salmon_versions=\$(sed -n 's/.*"salmon_version": *"\\([^"]*\\)".*/\\1/p' quant/*/aux_info/meta_info.json | sort -u)
    if [ "\$(echo "\$salmon_versions" | wc -l)" -ne 1 ] || [ -z "\$salmon_versions" ]; then
        echo "expected one Salmon version across samples, found: \$salmon_versions" >&2
        exit 1
    fi
    r_versions=\$(Rscript -e 'cat(sprintf("tximport_%s_DESeq2_%s", packageVersion("tximport"), packageVersion("DESeq2")))')
    echo "Salmon_\${salmon_versions}_\${r_versions}_unshrunk_effect" > deseq2/analysis_method.txt
    """

    stub:
    """
    mkdir deseq2 && touch deseq2/${params.test_level}_vs_${params.reference_level}_deseq2_all_genes.tsv
    echo "stub_analysis_method" > deseq2/analysis_method.txt
    """
}

process EXPORT_EVIDENCE {
    label 'aree'
    publishDir "${params.outdir}/evidence", mode: 'copy'

    input:
    path results
    path analysis_method
    path gff

    output:
    path "${params.study_id}_rnaseq.tsv", emit: processed
    path "${params.study_id}_mapping.tsv", emit: mapping
    path "${params.study_id}_annotations.tsv", emit: annotations

    script:
    def flags = params.quality_flags.tokenize(',').collect { "--quality-flag ${it.trim()}" }.join(' ')
    def comparison = params.sample_comparison ?: "${params.test_level}_vs_${params.reference_level}"
    """
    aree export-deseq2-evidence --results ${results} --gff ${gff} \\
        --processed ${params.study_id}_rnaseq.tsv --mapping ${params.study_id}_mapping.tsv \\
        --annotations ${params.study_id}_annotations.tsv --sample-comparison ${comparison} \\
        --reference-release '${params.reference_release}' --analysis-method "\$(cat ${analysis_method})" ${flags}
    """

    stub:
    """
    touch ${params.study_id}_rnaseq.tsv ${params.study_id}_mapping.tsv ${params.study_id}_annotations.tsv
    """
}

workflow {
    required('study_id')
    required('test_level')
    if (!params.salmon_index && !params.transcripts) {
        error "Provide --salmon_index (prebuilt) or --transcripts (optionally with --decoys) to build one"
    }

    reads = Channel
        .fromPath(required('samplesheet'), checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2) {
                error "Sample sheet rows need sample, fastq_1 and fastq_2: ${row}"
            }
            tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }
    design = file(required('design'), checkIfExists: true)
    tx2gene = file(required('tx2gene'), checkIfExists: true)
    gff = file(required('gff'), checkIfExists: true)

    if (!params.skip_fastqc) {
        FASTQC(reads)
    }
    index = params.salmon_index
        ? Channel.value(file(params.salmon_index, checkIfExists: true))
        : SALMON_INDEX(
            file(params.transcripts, checkIfExists: true),
            params.decoys ? file(params.decoys, checkIfExists: true) : file("${projectDir}/assets/NO_DECOYS")
        )
    quant = SALMON_QUANT(reads, index).collect()
    SALMON_QC(quant, design)
    DESEQ2(quant, design, tx2gene)
    EXPORT_EVIDENCE(DESEQ2.out.results, DESEQ2.out.analysis_method, gff)
}
