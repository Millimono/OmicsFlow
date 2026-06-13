/*
 * MODULE: FastQC
 * Performs quality control on raw FASTQ reads
 */

process FASTQC {

    tag "$meta.id"
    label 'process_medium'

    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir "${params.outdir}/qc/fastqc", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.zip') ? null : filename }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    path "versions.yml",             emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def mem_mb   = task.memory.toMega()
    """
    fastqc \\
        --threads ${task.cpus} \\
        --memory  ${mem_mb} \\
        --outdir  . \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_fastqc.html
    touch ${meta.id}_fastqc.zip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}
