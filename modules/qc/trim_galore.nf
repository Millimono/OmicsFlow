/*
 * MODULE: Trim Galore
 * Adapter trimming and quality filtering for Illumina reads
 */

process TRIM_GALORE {

    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'

    publishDir "${params.outdir}/trimmed", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.fq.gz'))     "reads/$filename"
            else if (filename.endsWith('_trimming_report.txt')) "logs/$filename"
            else null
        }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{_val_1,_val_2,_trimmed}.fq.gz"), emit: reads
    tuple val(meta), path("*_trimming_report.txt"),            emit: log
    path "versions.yml",                                       emit: versions

    script:
    def paired     = meta.single_end ? "" : "--paired"
    def min_len    = params.min_length    ?: 20
    def quality    = params.quality_cutoff ?: 20
    def cores      = task.cpus > 4 ? 4 : task.cpus   // Trim Galore max 4 cores
    """
    trim_galore \\
        ${paired} \\
        --quality    ${quality} \\
        --length     ${min_len} \\
        --cores      ${cores} \\
        --gzip \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(trim_galore --version | grep 'version' | sed 's/.*version //')
        cutadapt:   \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def suffix = meta.single_end ? "_trimmed.fq.gz" : "_val_1.fq.gz _val_2.fq.gz"
    """
    touch ${meta.id}${suffix}
    touch ${meta.id}_trimming_report.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: 0.6.7
        cutadapt:   3.4
    END_VERSIONS
    """
}
