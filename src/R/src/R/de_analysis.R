#' Perform differential expression analysis
#'
#' This function performs differential expression analysis using DESeq2.
#'
#' @param counts_file Path to the counts file
#' @param coldata_file Path to the coldata file
#' @param output_dir Directory to save results
#' @param lfc_thresh Log fold change threshold
#' @param pval_thresh P-value threshold
#' @param fdr False discovery rate
#' @param target Target variable name
#' @param reference_level Reference level for the condition
#' @param contrast_level Contrast level for the condition
#' @param covariates List of covariates
#' @export
de_analysis <- function(counts_file, coldata_file, output_dir, target,
                        reference_level, contrast_level, covariates,
                        lfc_thresh = 1, pval_thresh = 0.05,
                        fdr = 0.05) {
    suppressPackageStartupMessages({
        library(DESeq2)
        library(ggpubr)
        library(ggrepel)
        library(BiocParallel)
        library(tibble)
    })

    n_cores <- max(1, parallelly::availableCores() - 2)
    BiocParallel::register(MulticoreParam(n_cores))

    # Create output directory
    dir.create(output_dir, showWarnings = FALSE)

    mean_difference <- function(fit, name) {
        mean_diff_rna_seq_path <- file.path(output_dir, name)
        diff_expr <- data.frame(fit)
        options(ggrepel.max.overlaps = Inf)
        ggpubr::ggmaplot(diff_expr,
            main = expression("Group 1" %->% "   Group 2"),
            fdr = fdr, fc = 2, size = 0.4,
            palette = c("#B31B21", "#1465AC", "darkgray"),
            legend = "top", top = 0,
            font.legend = "bold", font.main = "bold",
            select.top.method = c("padj", "fc"),
            xlab = "Log2 mean expression", ylab = "Log2 fold change",
            ggtheme = ggplot2::theme_minimal()
        )
        ggsave(mean_diff_rna_seq_path, width = 15, height = 15, dpi = 300, units = "cm", bg = "white")
    }

    volcano_plot <- function(fit, name) {
        rna_volcano_path <- file.path(output_dir, name)
        de <- fit
        de$diff_expr <- "Not sig."
        de$diff_expr[de$log2FoldChange > lfc_thresh & de$pvalue < pval_thresh] <- "Up"
        de$diff_expr[de$log2FoldChange < -lfc_thresh & de$pvalue < pval_thresh] <- "Down"

        plot <- ggplot(data = de, aes(x = de$log2FoldChange, y = -log10(de$pvalue), col = diff_expr)) +
            geom_point() +
            theme_minimal() +
            scale_color_manual(values = c("#1465AC", "darkgray", "#B31B21")) +
            geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), col = "#B31B21") +
            geom_hline(yintercept = -log10(pval_thresh), col = "#B31B21")
        ggsave(rna_volcano_path, width = 15, height = 15, dpi = 300, units = "cm", bg = "white")
    }

    # Read counts
    if (grepl("\\.gz$", counts_file)) {
        counts <- read.table(gzfile(counts_file), header = TRUE, sep = ",", row.names = 1)
    } else {
        counts <- read.table(counts_file, header = TRUE, sep = ",", row.names = 1)
    }
    colnames(counts) <- gsub("\\.", "-", colnames(counts))

    # Read coldata
    coldata <- read.csv(coldata_file, header = TRUE, row.names = 1)
    colnames(coldata) <- gsub(paste0("^", target, "$"), "condition", colnames(coldata))

    # Convert covariates to factors, add covariates with >= 2 levels to the design formula
    for (covariate in covariates) {
        coldata[[covariate]] <- as.factor(coldata[[covariate]])
    }
    design_formula <- "~ "
    for (covariate in covariates) {
        if (length(levels(coldata[[covariate]])) > 1) {
            design_formula <- paste0(design_formula, covariate, " + ")
        }
    }
    design_formula <- paste(design_formula, "condition")
    design <- as.formula(design_formula)

    levels <- c(reference_level, contrast_level)
    coldata <- coldata[coldata$condition %in% levels, ]
    coldata$condition <- factor(coldata$condition, levels = levels)
    counts <- counts[, rownames(coldata)]

    # Filter low-count genes
    keep <- rowSums(counts > 5) >= (ncol(counts) / 2)
    counts <- counts[keep, ]

    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = design
    )
    dds <- DESeq(dds, parallel = TRUE)
    res <- results(
        dds,
        contrast = c("condition", contrast_level, reference_level),
        pAdjustMethod = "BH",
        alpha = fdr,
        parallel = TRUE
    )
    res <- as_tibble(res, rownames = "gene")

    cohort_name <- strsplit(basename(coldata_file), "\\.")[[1]][1]
    output_file <- file.path(output_dir, paste0("DE_results_", cohort_name, ".csv"))

    write.csv(as.data.frame(res), file = output_file, row.names = FALSE, quote = FALSE)

    volcano_plot(res, paste0("volcano_plot_", cohort_name, ".png"))
    mean_difference(res, paste0("mean_difference_", cohort_name, ".png"))

    cat("Saved results to", output_file, "\n")
}
