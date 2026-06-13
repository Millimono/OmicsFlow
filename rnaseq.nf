#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OmicsFlow — RNA-seq Pipeline
    Author  : Sory Millimono
    Version : 1.0.0
    License : MIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Usage:
        nextflow run rnaseq.nf \
            --input  data/test/samplesheet.csv \
            --genome GRCh38 \
            --outdir results/ \
            -profile docker

    Input samplesheet format (CSV):
        sample,fastq_1,fastq_2,strandedness
        sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,forward
        sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,reverse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ── PARAMETERS ──────────────────────────────────────────────────────────────
params {

    // Input / Output
    input         = null
    outdir        = "results"
    genome        = "GRCh38"

    // Reference files (auto-resolved if using iGenomes)
    fasta         = null
    gtf           = null
    star_index    = null
    salmon_index  = null

    // Trimming
    skip_trimming = false
    min_length    = 20
    quality_cutoff= 20

    // Alignment
    aligner       = "star"     // star | hisat2
    save_unaligned= false

    // Quantification
    quantifier    = "salmon"   // salmon | featurecounts

    // Differential expression
    skip_deseq2   = false
    min_counts    = 10
    fdr_cutoff    = 0.05
    lfc_threshold = 1.0

    // QC
    skip_fastqc   = false
    skip_multiqc  = false

    // Resources
    max_cpus      = 8
    max_memory    = "32.GB"
    max_time      = "24.h"
}

// ── INCLUDE MODULES ─────────────────────────────────────────────────────────
include { FASTQC         } from './modules/qc/fastqc'
include { TRIM_GALORE    } from './modules/qc/trim_galore'
include { STAR_ALIGN     } from './modules/alignment/star'
include { SALMON_QUANT   } from './modules/quantification/salmon'
include { MULTIQC        } from './modules/qc/multiqc'
include { SAMTOOLS_SORT  } from './modules/alignment/samtools'
include { DESEQ2_ANALYSIS} from './modules/quantification/deseq2'

// ── HELPER FUNCTIONS ────────────────────────────────────────────────────────

/*
 * Parse samplesheet CSV and return a channel of tuples:
 * [ meta, [ fastq_1, fastq_2 ] ]
 */
def parseSamplesheet(csv_file) {
    Channel
        .fromPath(csv_file)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id          : row.sample,
                strandedness: row.strandedness ?: 'unstranded',
                single_end  : !row.fastq_2 || row.fastq_2.isEmpty()
            ]
            def reads = meta.single_end
                ? [ file(row.fastq_1, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true) ]
            return [ meta, reads ]
        }
}

// ── MAIN WORKFLOW ────────────────────────────────────────────────────────────
workflow {

    // ── 0. Validate inputs ──────────────────────────────────────────────────
    if (!params.input) {
        error "Please provide a samplesheet: --input samplesheet.csv"
    }

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║           OmicsFlow — RNA-seq Pipeline               ║
    ╠══════════════════════════════════════════════════════╣
    ║  Input        : ${params.input}
    ║  Genome       : ${params.genome}
    ║  Aligner      : ${params.aligner}
    ║  Quantifier   : ${params.quantifier}
    ║  Output dir   : ${params.outdir}
    ║  Skip trimming: ${params.skip_trimming}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()

    // ── 1. Parse samplesheet ────────────────────────────────────────────────
    ch_reads = parseSamplesheet(params.input)

    // ── 2. FastQC on raw reads ──────────────────────────────────────────────
    ch_fastqc_raw = Channel.empty()
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_fastqc_raw = FASTQC.out.zip
    }

    // ── 3. Trimming ─────────────────────────────────────────────────────────
    ch_trimmed_reads = ch_reads
    ch_trim_log      = Channel.empty()

    if (!params.skip_trimming) {
        TRIM_GALORE(ch_reads)
        ch_trimmed_reads = TRIM_GALORE.out.reads
        ch_trim_log      = TRIM_GALORE.out.log
    }

    // ── 4. FastQC on trimmed reads ──────────────────────────────────────────
    ch_fastqc_trimmed = Channel.empty()
    if (!params.skip_fastqc && !params.skip_trimming) {
        FASTQC(ch_trimmed_reads.map { meta, reads -> [ meta + [trimmed: true], reads ] })
        ch_fastqc_trimmed = FASTQC.out.zip
    }

    // ── 5. Alignment ─────────────────────────────────────────────────────────
    ch_star_index = params.star_index
        ? Channel.value(file(params.star_index))
        : Channel.value([])

    STAR_ALIGN(ch_trimmed_reads, ch_star_index)
    ch_bam = STAR_ALIGN.out.bam

    // ── 6. Sort & index BAM ──────────────────────────────────────────────────
    SAMTOOLS_SORT(ch_bam)
    ch_bam_sorted = SAMTOOLS_SORT.out.bam

    // ── 7. Quantification ────────────────────────────────────────────────────
    ch_salmon_index = params.salmon_index
        ? Channel.value(file(params.salmon_index))
        : Channel.value([])

    SALMON_QUANT(ch_trimmed_reads, ch_salmon_index)
    ch_counts = SALMON_QUANT.out.results

    // ── 8. Differential expression ───────────────────────────────────────────
    if (!params.skip_deseq2) {
        ch_counts_collected = ch_counts.collect()
        DESEQ2_ANALYSIS(ch_counts_collected)
    }

    // ── 9. MultiQC ───────────────────────────────────────────────────────────
    if (!params.skip_multiqc) {
        ch_multiqc_files = Channel.empty()
            .mix(ch_fastqc_raw)
            .mix(ch_fastqc_trimmed)
            .mix(ch_trim_log)
            .mix(STAR_ALIGN.out.log)
            .mix(SALMON_QUANT.out.log)
            .collect()

        MULTIQC(ch_multiqc_files)
    }

    // ── 10. Summary ──────────────────────────────────────────────────────────
    workflow.onComplete {
        log.info """
        ╔══════════════════════════════════════════════════════╗
        ║              Pipeline completed!                     ║
        ╠══════════════════════════════════════════════════════╣
        ║  Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        ║  Duration  : ${workflow.duration}
        ║  Results   : ${params.outdir}/
        ║  Report    : ${params.outdir}/pipeline_report.html
        ╚══════════════════════════════════════════════════════╝
        """.stripIndent()
    }
}

// ── WORKFLOW: QC ONLY (subworkflow) ─────────────────────────────────────────
workflow QC_ONLY {
    take: ch_reads
    main:
        FASTQC(ch_reads)
        if (!params.skip_trimming) {
            TRIM_GALORE(ch_reads)
            FASTQC(TRIM_GALORE.out.reads.map { m, r -> [ m + [trimmed: true], r ] })
        }
        MULTIQC(FASTQC.out.zip.collect())
    emit:
        report = MULTIQC.out.report
}
