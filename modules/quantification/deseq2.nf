/*
 * MODULE: DESeq2
 * Differential expression analysis from Salmon counts
 */

process DESEQ2_ANALYSIS {

    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_0'

    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path salmon_dirs   // collected list of Salmon output dirs

    output:
    path "deseq2_results.csv",   emit: results
    path "normalized_counts.csv",emit: counts
    path "volcano_plot.pdf",     emit: volcano
    path "heatmap_top50.pdf",    emit: heatmap
    path "pca_plot.pdf",         emit: pca
    path "versions.yml",         emit: versions

    script:
    def fdr = params.fdr_cutoff    ?: 0.05
    def lfc = params.lfc_threshold ?: 1.0
    def min = params.min_counts    ?: 10
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(DESeq2)
        library(tximport)
        library(ggplot2)
        library(pheatmap)
        library(RColorBrewer)
    })

    # ── 1. Import Salmon counts via tximport ──────────────────────────────
    salmon_files <- list.files(
        path    = ".",
        pattern = "quant.sf",
        recursive = TRUE,
        full.names = TRUE
    )
    names(salmon_files) <- basename(dirname(salmon_files))

    txi <- tximport(
        files  = salmon_files,
        type   = "salmon",
        txOut  = FALSE   # summarize to gene level
    )

    # ── 2. Build colData from sample names ────────────────────────────────
    sample_names <- names(salmon_files)
    # Simple condition parsing: sample names ending in _ctrl or _treat
    condition <- ifelse(grepl("ctrl|control|untreated", sample_names, ignore.case = TRUE),
                        "control", "treatment")
    coldata <- data.frame(
        row.names = sample_names,
        condition = factor(condition, levels = c("control", "treatment"))
    )
    message("Samples detected: ", paste(sample_names, collapse=", "))
    message("Conditions: ", paste(condition, collapse=", "))

    # ── 3. DESeq2 ─────────────────────────────────────────────────────────
    dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

    # Filter low-count genes
    keep <- rowSums(counts(dds)) >= ${min}
    dds  <- dds[keep, ]
    message("Genes after filtering: ", nrow(dds))

    dds  <- DESeq(dds)
    res  <- results(dds,
                    contrast  = c("condition", "treatment", "control"),
                    alpha     = ${fdr},
                    lfcThreshold = ${lfc})
    res  <- lfcShrink(dds, coef = "condition_treatment_vs_control",
                       type = "apeglm")

    # ── 4. Export results ─────────────────────────────────────────────────
    res_df <- as.data.frame(res)
    res_df\$gene <- rownames(res_df)
    res_df <- res_df[order(res_df\$padj, na.last = TRUE), ]
    write.csv(res_df, "deseq2_results.csv", row.names = FALSE)

    norm_counts <- counts(dds, normalized = TRUE)
    write.csv(norm_counts, "normalized_counts.csv")

    message("Significant genes (padj < ${fdr}, |LFC| > ${lfc}): ",
            sum(res_df\$padj < ${fdr} & abs(res_df\$log2FoldChange) > ${lfc}, na.rm = TRUE))

    # ── 5. Volcano plot ───────────────────────────────────────────────────
    res_df\$sig <- ifelse(
        !is.na(res_df\$padj) & res_df\$padj < ${fdr} & abs(res_df\$log2FoldChange) > ${lfc},
        ifelse(res_df\$log2FoldChange > 0, "UP", "DOWN"), "NS"
    )
    counts_sig <- table(res_df\$sig)

    p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_manual(
            values = c("UP" = "#E74C3C", "DOWN" = "#3498DB", "NS" = "#BDC3C7"),
            labels = c(
                paste0("UP (n=",   counts_sig["UP"]   %||% 0, ")"),
                paste0("DOWN (n=", counts_sig["DOWN"] %||% 0, ")"),
                paste0("NS (n=",   counts_sig["NS"]   %||% 0, ")")
            )
        ) +
        geom_hline(yintercept = -log10(${fdr}), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-${lfc}, ${lfc}), linetype = "dashed", color = "grey50") +
        labs(
            title    = "Differential Expression — Volcano Plot",
            subtitle = paste0("FDR < ${fdr}  |  |LFC| > ${lfc}"),
            x        = "Log2 Fold Change",
            y        = "-Log10 adjusted p-value",
            color    = "Regulation"
        ) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("volcano_plot.pdf", p_volcano, width = 8, height = 6)

    # ── 6. Heatmap top 50 DE genes ────────────────────────────────────────
    sig_genes <- rownames(res_df)[!is.na(res_df\$padj) & res_df\$padj < ${fdr}]
    top_genes <- head(sig_genes, 50)

    if (length(top_genes) >= 2) {
        mat <- assay(vst(dds, blind = FALSE))[top_genes, ]
        mat <- mat - rowMeans(mat)
        pdf("heatmap_top50.pdf", width = 10, height = 12)
        pheatmap(
            mat,
            color            = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
            annotation_col   = coldata,
            show_rownames    = length(top_genes) <= 50,
            cluster_rows     = TRUE,
            cluster_cols     = TRUE,
            fontsize         = 8,
            main             = "Top DE genes (VST normalized, centered)"
        )
        dev.off()
    } else {
        pdf("heatmap_top50.pdf"); plot.new()
        text(0.5, 0.5, "Not enough significant genes for heatmap", cex = 1.2)
        dev.off()
    }

    # ── 7. PCA plot ───────────────────────────────────────────────────────
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"))

    p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -0.8, size = 3) +
        xlab(paste0("PC1: ", pct_var[1], "% variance")) +
        ylab(paste0("PC2: ", pct_var[2], "% variance")) +
        labs(title = "PCA — Sample clustering") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave("pca_plot.pdf", p_pca, width = 7, height = 6)

    message("DESeq2 analysis complete.")

    # versions
    writeLines(
        c('"${task.process}":', paste0('    deseq2: ', packageVersion("DESeq2"))),
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
    END_VERSIONS
    """
}
