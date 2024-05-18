#' Perform variance stabilizing transformation (VST)
#'
#' This function performs VST on RNA-seq count data.
#'
#' @param counts_file Path to the counts file
#' @param coldata_file Path to the coldata file
#' @param output_dir Directory to save VST results
#' @param design_formula_cols Coldata fields to include in design formula
#' @export
vst_transform <- function(counts_file, coldata_file, output_dir, design_formula_cols = NULL) {
    suppressPackageStartupMessages(library(DESeq2))

    # Create output directory
    dir.create(output_dir, showWarnings = FALSE)

    # Read counts table
    if (grepl("\\.gz$", counts_file)) {
        counts <- read.table(gzfile(counts_file), header = TRUE, sep = ",", row.names = 1)
    } else {
        counts <- read.table(counts_file, header = TRUE, sep = ",", row.names = 1)
    }
    colnames(counts) <- gsub("\\.", "-", colnames(counts))

    # Read coldata
    coldata <- read.csv(coldata_file, header = TRUE, row.names = 1)

    if (!is.null(design_formula_cols)) {
        # Convert covariates to factors
        for (covariate in design_formula_cols) {
            coldata[[covariate]] <- as.factor(coldata[[covariate]])
        }

        # Add covariates with >= 2 levels to the design formula
        design_formula <- "~ "
        for (covariate in design_formula_cols) {
            if (length(levels(coldata[[covariate]])) > 1) {
                design_formula <- paste0(design_formula, covariate, " + ")
            }
        }
        if (design_formula == "~ ") {
            design_formula <- "~ 1"
        } else {
            # Remove trailing " + " from design formula
            design_formula <- substr(design_formula, 1, nchar(design_formula) - 3)
        }
    } else {
        design_formula <- "~ 1"
    }
    cat("Design formula:", design_formula, "\n")

    dds <- DESeqDataSetFromMatrix(
        countData = counts, colData = coldata,
        design = as.formula(design_formula)
    )

    # Perform variance stabilizing transformation
    vst_data <- vst(dds)
    transformed_data <- assay(vst_data)

    # Prepare transformed data for output
    transformed_data <- data.frame(Ensembl_gene_id = rownames(transformed_data), transformed_data)
    colnames(transformed_data) <- gsub("\\.", "-", colnames(transformed_data))

    # Write output
    output_path <- file.path(output_dir, basename(counts_file))
    write.table(transformed_data, gzfile(output_path), sep = ",", quote = FALSE, row.names = FALSE)

    cat("VST transformed data written to:", output_path, "\n")
}
