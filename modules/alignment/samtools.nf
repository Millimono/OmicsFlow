/*
 * MODULE: Samtools sort + index
 */

process SAMTOOLS_SORT {

    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'

    publishDir "${params.outdir}/aligned/bam", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"),     emit: bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: bai
    tuple val(meta), path("*.flagstat"),        emit: flagstat
    path "versions.yml",                        emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def mem_per_thread = "${(task.memory.toGiga() / task.cpus).toInteger()}G"
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${mem_per_thread} \\
        -o ${prefix}.sorted.bam \\
        ${bam}

    samtools index ${prefix}.sorted.bam

    samtools flagstat ${prefix}.sorted.bam > ${prefix}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.flagstat
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.18
    END_VERSIONS
    """
}
