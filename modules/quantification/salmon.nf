/*
 * MODULE: Salmon
 * Transcript-level quantification (alignment-free, fast)
 */

process SALMON_QUANT {

    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/salmon:1.10.2--hecfa306_0'

    publishDir "${params.outdir}/counts/salmon", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${meta.id}"),             emit: results
    tuple val(meta), path("${meta.id}/aux_info/"),   emit: aux
    tuple val(meta), path("${meta.id}/logs/"),        emit: log
    path "versions.yml",                              emit: versions

    script:
    def lib_type  = meta.strandedness == 'forward'  ? 'SF' :
                    meta.strandedness == 'reverse'   ? 'SR' :
                    meta.single_end                  ? 'U'  : 'IU'
    def reads_arg = meta.single_end
                    ? "-r ${reads}"
                    : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    salmon quant \\
        --index      ${index} \\
        --libType    ${lib_type} \\
        ${reads_arg} \\
        --threads    ${task.cpus} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        --output     ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/aux_info ${meta.id}/logs
    touch ${meta.id}/quant.sf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: 1.10.2
    END_VERSIONS
    """
}
