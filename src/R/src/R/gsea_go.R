#' Perform GSEA GO analysis
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using Gene Ontology (GO) terms.
#'
#' @param de_file Path to the DE results file
#' @param gsea_directory Directory to save GSEA results
#' @export
gsea_go <- function(de_file, gsea_directory) {
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(tibble)
        library(tidyr)
        library(BiocParallel)
        library(AnnotationDbi)
        library(org.Hs.eg.db)
        library(clusterProfiler)
        library(rrvgo)
        library(ggplot2)
        library(plotly)
        library(enrichplot)
    })

    n_cores <- max(1, parallelly::availableCores() - 2)
    BiocParallel::register(MulticoreParam(n_cores))

    dir.create(gsea_directory, showWarnings = FALSE)

    set.seed(0)

    experiment_name <- gsub("DE_results_", "", basename(de_file))
    experiment_name <- gsub(".csv", "", experiment_name)
    cat("Processing:", experiment_name, "\n")
    res <- fread(de_file)

    ens2entrez <- AnnotationDbi::select(org.Hs.eg.db,
        keys = res$gene,
        columns = "ENTREZID",
        keytype = "ENSEMBL"
    )

    res <- inner_join(res, as_tibble(ens2entrez), by = c("gene" = "ENSEMBL"))

    res <- filter(res, !is.na(ENTREZID) & !is.na(stat))
    res <- res[!duplicated(res$ENTREZID), ]

    values <- res$stat
    names(values) <- res$ENTREZID
    values <- sort(values, decreasing = TRUE)

    # Using GO biological process ontology
    gseaGO_BP <- gseGO(geneList = values,
                        ont = "BP",
                        keyType = "ENTREZID",
                        OrgDb = org.Hs.eg.db,
                        pvalueCutoff = 0.05,
                        verbose = TRUE,
                        eps = 0)
    gseaGO_BP <- clusterProfiler::simplify(gseaGO_BP)

    go_analysis <- gseaGO_BP@result

    # Top 35 terms
    go_analysis_top_terms <- go_analysis[1:35, ]

    # Using rrvgo here, vignette:
    # https://bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html#getting-started
    # rrvgo helps deal with the redundancy between GO terms and for visualization

    simMatrix <- calculateSimMatrix(
        go_analysis_top_terms$ID,
        orgdb = "org.Hs.eg.db",
        ont = "BP",
        method = "Rel"
    )

    scores <- setNames(-log10(go_analysis_top_terms$qvalue), go_analysis_top_terms$ID)
    reducedTerms <- reduceSimMatrix(
        simMatrix,
        scores,
        threshold = 0.6,
        orgdb = "org.Hs.eg.db"
    )

    # Save heatmap, treemap, and scatterplot
    png(file.path(gsea_directory, paste0("heatmap_", experiment_name, ".png")),
        width = 1200, height = 1200
    )
    heatmapPlot(simMatrix,
        reducedTerms,
        annotateParent = TRUE,
        annotationLabel = "parentTerm",
        fontsize = 10
    )
    dev.off()

    png(file.path(gsea_directory, paste0("treemap_", experiment_name, ".png")))
    treemapPlot(reducedTerms)
    dev.off()

    top_terms <- unique(reducedTerms$parent)
    go_analysis_top_terms <- go_analysis_top_terms[go_analysis_top_terms$ID %in% top_terms, ]

    # convert Entrez IDs (list of genes separted by "/") to symbol
    go_analysis_top_terms$core_enrichment_symbols <- sapply(strsplit(go_analysis_top_terms$core_enrichment, "/"), function(x) {
        symbols <- AnnotationDbi::select(org.Hs.eg.db,
            keys = x,
            columns = "SYMBOL",
            keytype = "ENTREZID"
        )
        paste(symbols$SYMBOL, collapse = "/")
    })

    # convert Entrez IDs (list of genes separated by "/") to Ensembl IDs
    go_analysis_top_terms$core_enrichment_ensembl <- sapply(strsplit(go_analysis_top_terms$core_enrichment, "/"), function(x) {
        ensembl <- AnnotationDbi::select(org.Hs.eg.db,
            keys = x,
            columns = "ENSEMBL",
            keytype = "ENTREZID"
        )
        paste(ensembl$ENSEMBL, collapse = "/")
    })

    write.csv(go_analysis, file = file.path(gsea_directory, paste0("gsea_go_bp_", experiment_name, ".csv")), row.names = FALSE)
    write.csv(go_analysis_top_terms, file = file.path(gsea_directory, paste0("gsea_go_bp_top_terms_", experiment_name, ".csv")), row.names = FALSE)

    # Generate interactive Plotly scatter plot
    simMatrix <- simMatrix[rownames(simMatrix) %in% top_terms, ]
    reducedTerms <- reducedTerms[reducedTerms$go %in% top_terms, ]
    scatter <- scatterPlot(simMatrix, reducedTerms, size = "size")
    ggsave(file.path(gsea_directory, paste0("scatter_", experiment_name, ".png")), plot = scatter)
    p <- ggplotly(scatter)
    htmlwidgets::saveWidget(p, file.path(gsea_directory, paste0("scatter_", experiment_name, ".html")))
}
