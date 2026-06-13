/*
 * MODULE: STAR
 * Splice-aware alignment of RNA-seq reads to a reference genome
 */

process STAR_ALIGN {

    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/star:2.7.11a--h0033a41_0'

    publishDir "${params.outdir}/aligned", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.bam'))       "bam/$filename"
            else if (filename.endsWith('Log.final.out')) "logs/$filename"
            else if (filename.endsWith('SJ.out.tab'))    "junctions/$filename"
            else null
        }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*.Log.final.out"),                 emit: log
    tuple val(meta), path("*.SJ.out.tab"),                    emit: junctions
    tuple val(meta), path("*.Unmapped*"),     optional: true, emit: unmapped
    path "versions.yml",                                      emit: versions

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def strand_field = meta.strandedness == 'forward'  ? '1' :
                       meta.strandedness == 'reverse'  ? '2' : '0'
    def reads_arg    = meta.single_end
                        ? reads
                        : "${reads[0]} ${reads[1]}"
    def unmapped_arg = params.save_unaligned
                        ? "--outReadsUnmapped Fastx"
                        : ""
    """
    STAR \\
        --runMode            alignReads \\
        --runThreadN         ${task.cpus} \\
        --genomeDir          ${index} \\
        --readFilesIn        ${reads_arg} \\
        --readFilesCommand   zcat \\
        --outSAMtype         BAM SortedByCoordinate \\
        --outSAMattributes   NH HI AS NM MD \\
        --outSAMstrandField  intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonical \\
        --outFileNamePrefix  ${prefix}. \\
        --sjdbScore          1 \\
        --quantMode          TranscriptomeSAM \\
        --twopassMode        Basic \\
        ${unmapped_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.SJ.out.tab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11a
    END_VERSIONS
    """
}
