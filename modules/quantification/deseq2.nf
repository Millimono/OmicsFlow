/*
 * MODULE: DESeq2
 * Differential expression analysis from Salmon counts
 * Fix 1: tx2gene support for tximport gene-level summarization
 * Fix 2: conditions read from samplesheet (not inferred from names)
 */

process DESEQ2_ANALYSIS {

    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_0'

    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path salmon_dirs       // collected list of Salmon output dirs
    path samplesheet       // CSV with sample,condition columns
    path tx2gene           // transcript-to-gene mapping (optional)

    output:
    path "deseq2_results.csv",    emit: results
    path "normalized_counts.csv", emit: counts
    path "volcano_plot.pdf",      emit: volcano
    path "heatmap_top50.pdf",     emit: heatmap
    path "pca_plot.pdf",          emit: pca
    path "versions.yml",          emit: versions

    script:
    def fdr    = params.fdr_cutoff    ?: 0.05
    def lfc    = params.lfc_threshold ?: 1.0
    def min    = params.min_counts    ?: 10
    def tx2g   = tx2gene.name != 'NO_FILE' ? "TRUE" : "FALSE"
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(DESeq2)
        library(tximport)
        library(ggplot2)
        library(pheatmap)
        library(RColorBrewer)
    })

    # ── 1. Read samplesheet for conditions ──────────────────────────────────
    # Expected columns: sample, condition (at minimum)
    meta <- read.csv("${samplesheet}", stringsAsFactors = FALSE)

    # Validate required columns
    required_cols <- c("sample", "condition")
    missing <- setdiff(required_cols, colnames(meta))
    if (length(missing) > 0) {
        stop("Samplesheet is missing required columns: ",
             paste(missing, collapse = ", "),
             "\\nExpected: sample,fastq_1,fastq_2,strandedness,condition")
    }

    rownames(meta) <- meta\$sample
    message("Samples: ", paste(meta\$sample, collapse = ", "))
    message("Conditions: ", paste(meta\$condition, collapse = ", "))

    # Validate at least 2 conditions
    if (length(unique(meta\$condition)) < 2) {
        stop("At least 2 distinct conditions required. Found: ",
             paste(unique(meta\$condition), collapse = ", "))
    }

    # ── 2. Locate Salmon quant.sf files ─────────────────────────────────────
    salmon_files <- list.files(
        path      = ".",
        pattern   = "quant.sf",
        recursive = TRUE,
        full.names = TRUE
    )
    names(salmon_files) <- basename(dirname(salmon_files))

    # Reorder to match samplesheet
    common <- intersect(meta\$sample, names(salmon_files))
    if (length(common) == 0) {
        stop("No matching samples between samplesheet and Salmon outputs.\\n",
             "Samplesheet samples: ", paste(meta\$sample, collapse=", "), "\\n",
             "Salmon dirs found:   ", paste(names(salmon_files), collapse=", "))
    }
    salmon_files <- salmon_files[common]
    meta         <- meta[common, , drop = FALSE]

    # ── 3. tximport — with or without tx2gene ───────────────────────────────
    use_tx2gene <- ${tx2g}

    if (use_tx2gene) {
        message("Using tx2gene file for gene-level summarization")
        tx2gene <- read.csv("${tx2gene}", header = TRUE,
                            col.names = c("tx_id", "gene_id"))
        txi <- tximport(
            files    = salmon_files,
            type     = "salmon",
            tx2gene  = tx2gene,
            txOut    = FALSE
        )
    } else {
        message("No tx2gene provided — using transcript-level counts (txOut=TRUE)")
        message("For gene-level analysis, provide --tx2gene tx2gene.csv")
        txi <- tximport(
            files = salmon_files,
            type  = "salmon",
            txOut = TRUE
        )
    }

    # ── 4. DESeq2 ────────────────────────────────────────────────────────────
    coldata <- data.frame(
        row.names = common,
        condition = factor(meta\$condition)
    )
    # Set reference level (alphabetically first, or "control" if present)
    ref_level <- ifelse("control" %in% levels(coldata\$condition),
                        "control", levels(coldata\$condition)[1])
    coldata\$condition <- relevel(coldata\$condition, ref = ref_level)
    message("Reference condition: ", ref_level)

    dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

    # Filter low-count genes
    keep <- rowSums(counts(dds)) >= ${min}
    dds  <- dds[keep, ]
    message("Genes after filtering (>= ${min} counts): ", nrow(dds))

    dds <- DESeq(dds)

    # Get contrast: treatment vs reference
    contrast_levels <- levels(coldata\$condition)
    treat_level     <- contrast_levels[contrast_levels != ref_level][1]
    message("Contrast: ", treat_level, " vs ", ref_level)

    res <- results(dds,
                   contrast     = c("condition", treat_level, ref_level),
                   alpha        = ${fdr},
                   lfcThreshold = ${lfc})

    # LFC shrinkage with apeglm
    coef_name <- paste0("condition_", treat_level, "_vs_", ref_level)
    res <- tryCatch(
        lfcShrink(dds, coef = coef_name, type = "apeglm"),
        error = function(e) {
            message("lfcShrink failed (", e\$message, "), using raw LFC")
            res
        }
    )

    # ── 5. Export results ────────────────────────────────────────────────────
    res_df        <- as.data.frame(res)
    res_df\$gene  <- rownames(res_df)
    res_df        <- res_df[order(res_df\$padj, na.last = TRUE), ]
    write.csv(res_df, "deseq2_results.csv", row.names = FALSE)

    norm_counts <- counts(dds, normalized = TRUE)
    write.csv(norm_counts, "normalized_counts.csv")

    n_sig <- sum(res_df\$padj < ${fdr} & abs(res_df\$log2FoldChange) > ${lfc},
                 na.rm = TRUE)
    message("Significant DE genes (padj < ${fdr}, |LFC| > ${lfc}): ", n_sig)

    # ── 6. Volcano plot ──────────────────────────────────────────────────────
    res_df\$sig <- ifelse(
        !is.na(res_df\$padj) &
        res_df\$padj < ${fdr} &
        abs(res_df\$log2FoldChange) > ${lfc},
        ifelse(res_df\$log2FoldChange > 0, "UP", "DOWN"), "NS"
    )
    counts_sig <- table(factor(res_df\$sig, levels = c("UP","DOWN","NS")))

    p_volcano <- ggplot(res_df,
                        aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_manual(
            values = c("UP"="#E74C3C","DOWN"="#3498DB","NS"="#BDC3C7"),
            labels = c(paste0("UP (n=",   counts_sig["UP"],   ")"),
                       paste0("DOWN (n=", counts_sig["DOWN"], ")"),
                       paste0("NS (n=",   counts_sig["NS"],   ")"))
        ) +
        geom_hline(yintercept = -log10(${fdr}),
                   linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-${lfc}, ${lfc}),
                   linetype = "dashed", color = "grey50") +
        labs(title    = paste0("Differential Expression — ", treat_level, " vs ", ref_level),
             subtitle = paste0("FDR < ${fdr}  |  |LFC| > ${lfc}"),
             x = "Log2 Fold Change",
             y = "-Log10 adjusted p-value",
             color = "Regulation") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("volcano_plot.pdf", p_volcano, width = 8, height = 6)

    # ── 7. Heatmap top 50 DE genes ───────────────────────────────────────────
    sig_genes <- rownames(res_df)[!is.na(res_df\$padj) &
                                   res_df\$padj < ${fdr}]
    top_genes <- head(sig_genes, 50)

    if (length(top_genes) >= 2) {
        mat <- assay(vst(dds, blind = FALSE))[top_genes, ]
        mat <- mat - rowMeans(mat)
        annotation_col <- data.frame(
            condition = coldata\$condition,
            row.names = rownames(coldata)
        )
        pdf("heatmap_top50.pdf", width = 10, height = 12)
        pheatmap(
            mat,
            color          = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
            annotation_col = annotation_col,
            show_rownames  = length(top_genes) <= 50,
            cluster_rows   = TRUE,
            cluster_cols   = TRUE,
            fontsize       = 8,
            main           = "Top DE genes (VST normalized, centered)"
        )
        dev.off()
    } else {
        pdf("heatmap_top50.pdf"); plot.new()
        text(0.5, 0.5, "Not enough significant genes for heatmap", cex = 1.2)
        dev.off()
    }

    # ── 8. PCA plot ───────────────────────────────────────────────────────────
    vsd      <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 3) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        labs(title = paste0("PCA — Sample clustering (",
                            treat_level, " vs ", ref_level, ")")) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("pca_plot.pdf", p_pca, width = 7, height = 6)

    message("DESeq2 analysis complete.")

    # versions
    writeLines(
        c('"${task.process}":',
          paste0('    deseq2: ', packageVersion("DESeq2")),
          paste0('    tximport: ', packageVersion("tximport"))),
        "versions.yml"
    )
    """

    stub:
    """
    touch deseq2_results.csv normalized_counts.csv
    touch volcano_plot.pdf heatmap_top50.pdf pca_plot.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: 1.42.0
        tximport: 1.30.0
    END_VERSIONS
    """
}
