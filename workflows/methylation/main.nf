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

process BISMARK_GENOME_PREP {
    label 'bio'

    input:
    path genome_fasta

    output:
    path "bismark_index"

    script:
    """
    mkdir -p bismark_index
    cp ${genome_fasta} bismark_index/
    bismark_genome_preparation --bowtie2 bismark_index
    """

    stub:
    """
    mkdir -p bismark_index/Bisulfite_Genome
    """
}

process BISMARK_ALIGN_EXTRACT {
    tag "${sample}"
    label 'bio'
    publishDir "${params.outdir}/bismark", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path index_dir

    output:
    path "${sample}.bismark.cov.gz", emit: cov

    script:
    def read_arg = reads.size() > 1 ? "-1 ${reads[0]} -2 ${reads[1]}" : "${reads[0]}"
    """
    bismark --genome ${index_dir} ${read_arg} -o . --basename ${sample} --parallel ${task.cpus}
    bismark_methylation_extractor --bedGraph --gzip --counts --cytosine_report --no_header *.bam
    mv *.bismark.cov.gz ${sample}.bismark.cov.gz
    """

    stub:
    """
    echo -e "NC_TEST\\t100\\t100\\t50.0\\t5\\t5" | gzip -c > ${sample}.bismark.cov.gz
    """
}

process MERGE_METHYLATION_CALLS {
    label 'aree'
    publishDir "${params.outdir}/methylation_calls", mode: 'copy'

    input:
    path cov_files

    output:
    path "merged_methylation_coverage.tsv", emit: merged

    script:
    """
    python3 - <<'EOF'
    import gzip, glob, csv
    files = sorted(glob.glob("*.bismark.cov.gz"))
    samples = [f.replace(".bismark.cov.gz", "") for f in files]
    sites = {}
    for s, fn in zip(samples, files):
        with gzip.open(fn, "rt") as h:
            for line in h:
                parts = line.rstrip("\\n").split("\\t")
                if len(parts) < 6:
                    continue
                fid = f"{parts[0]}:{parts[1]}"
                sites.setdefault(fid, {})[s] = (int(parts[4]), int(parts[5]))
    fieldnames = ["feature_id"]
    for s in samples:
        fieldnames.extend([f"{s}_M", f"{s}_U"])
    with open("merged_methylation_coverage.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\\t")
        w.writeheader()
        for fid, smap in sorted(sites.items()):
            row = {"feature_id": fid}
            for s in samples:
                m, u = smap.get(s, (0, 0))
                row[f"{s}_M"] = m
                row[f"{s}_U"] = u
            w.writerow(row)
    EOF
    """

    stub:
    """
    echo -e "feature_id\\tcontrol_A_M\\tcontrol_A_U\\theat_A_M\\theat_A_U" > merged_methylation_coverage.tsv
    echo -e "NC_TEST:100\\t5\\t5\\t8\\t2" >> merged_methylation_coverage.tsv
    """
}

process LIMMA_DM {
    label 'bio'
    publishDir "${params.outdir}/differential", mode: 'copy'

    input:
    path coverage_table
    path design

    output:
    path "dm/*_differential_methylation.tsv", emit: results
    path "dm/analysis_method.txt", emit: analysis_method

    script:
    def covArg = params.covariate ? "--covariate=${params.covariate}" : ""
    """
    mkdir -p dm
    Rscript ${projectDir}/scripts/run_methylation_limma.R \\
        ${coverage_table} ${design} dm \\
        --reference=${params.reference_level} \\
        --test=${params.test_level} \\
        --feature-type=${params.feature_type} \\
        --min-coverage=${params.min_coverage} \\
        --min-samples=${params.min_samples} \\
        ${covArg}
    """
}

process EXPORT_METHYLATION_EVIDENCE {
    label 'aree'
    publishDir "${params.outdir}/evidence", mode: 'copy'

    input:
    path results
    path analysis_method

    output:
    path "${params.study_id}_methylation.tsv", emit: processed
    path "${params.study_id}_mapping.tsv", emit: mapping
    path "${params.study_id}_annotations.tsv", emit: annotations

    script:
    def flags = params.quality_flags.tokenize(',').collect { "--quality-flag ${shQuote(it.trim())}" }.join(' ')
    def comparison = params.sample_comparison ?: "${params.test_level}_vs_${params.reference_level}"
    """
    aree export-assay-evidence \\
        --results ${results} \\
        --processed ${params.study_id}_methylation.tsv \\
        --mapping ${params.study_id}_mapping.tsv \\
        --annotations ${params.study_id}_annotations.tsv \\
        --sample-comparison ${shQuote(comparison)} \\
        --feature-type ${shQuote(params.feature_type)} \\
        --effect-size-type m_value_diff \\
        --direction-mode methylation \\
        --reference-release ${shQuote(params.reference_release)} \\
        --analysis-method "\$(cat ${analysis_method})" \\
        ${flags}
    """
}

workflow {
    required('study_id')
    required('test_level')
    design = file(required('design'), checkIfExists: true)

    if (params.coverage_table) {
        cov = file(params.coverage_table, checkIfExists: true)
    } else if (params.samplesheet) {
        reads = Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def files = [file(row.fastq_1, checkIfExists: true)]
                if (row.fastq_2 && !row.fastq_2.toString().trim().isEmpty()) {
                    files << file(row.fastq_2, checkIfExists: true)
                }
                tuple(row.sample, files)
            }
        idx = params.bismark_index
            ? Channel.value(file(params.bismark_index, checkIfExists: true))
            : BISMARK_GENOME_PREP(file(required('genome_fasta'), checkIfExists: true))
        cov_calls = BISMARK_ALIGN_EXTRACT(reads, idx).cov.collect()
        cov = MERGE_METHYLATION_CALLS(cov_calls).merged
    } else {
        error "Provide --coverage_table (CpG/DMR matrix) or --samplesheet with Bismark genome"
    }

    LIMMA_DM(cov, design)
    EXPORT_METHYLATION_EVIDENCE(LIMMA_DM.out.results, LIMMA_DM.out.analysis_method)
    if (params.validate_evidence) {
        VALIDATE_EVIDENCE(
            EXPORT_METHYLATION_EVIDENCE.out.processed,
            EXPORT_METHYLATION_EVIDENCE.out.mapping,
            EXPORT_METHYLATION_EVIDENCE.out.annotations,
        )
    }
}
