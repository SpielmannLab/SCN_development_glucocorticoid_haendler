#!/usr/bin/env Rscript
# Call from conda-environment MegaBundle
" Make heatmap plots based on the GSEA analysis results

Usage: gsea_heatmap_per_pathway.R --file_sc_gsea_metadata=<file> --metadata_to_plot=<value> --heatmap_column=<value> --heatmap_row=<value>


Options:
  -h --help			Show this screen.
  --file_sc_gsea_metadata=<file>		RDS file containing sc_obj metadata and ssGSEA scores
  --metadata_to_plot=<value>  Use <ESCAPE_>. This has been appendend to all ssGSEA enrichment score colnames
  --heatmap_column=<value>  Which metadata column of Seurat object to use as columns of heatmap
  --heatmap_row=<value>     Which metadata column of Seurat object to use as rows of heatmap

" -> doc

# --- load libraries
suppressMessages(library(docopt))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressPackageStartupMessages(library(pheatmap))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

list2env(x = arguments, envir = environment())

message("file_sc_gsea_metadata: ", file_sc_gsea_metadata)
message("metadata_to_plot: ", metadata_to_plot)
message("heatmap_column: ", heatmap_column)
message("heatmap_row: ", heatmap_row)

# --- Read files
message("Loading Data")
sc_gsea_metadata <- readRDS(file_sc_gsea_metadata)


make_heatmap_per_pathway <- function(sc_gsea_metadata, heatmap_column, heatmap_row, pathway) {
    # get median gsea scores segregated by annotations for each of the Gene sets
    medians_all <- sc_gsea_metadata %>%
        select(all_of(c(heatmap_column, heatmap_row, pathway))) %>%
        group_by(across(all_of(c(heatmap_column, heatmap_row)))) %>%
        dplyr::summarise(across(all_of(pathway), median), .groups = "drop") %>%
        pivot_wider(names_from = all_of(heatmap_column), values_from = all_of(pathway)) %>%
        tibble::column_to_rownames(var = heatmap_row) %>%
        as.matrix() # convert to matrix for pheatmap

    pathway <- pathway %>%
        gsub(pattern = metadata_to_plot, replacement = "")

    filename = paste0("Heatmap_",
        gsub(basename(file_sc_gsea_metadata),
            pattern=".rds",
            replacement=paste0("_", heatmap_column, "_", heatmap_row, "_", pathway, ".pdf")))
    pheatmap(medians_all,
        cellheight = 20,
        cellwidth = 20,
        filename = filename,
        cluster_rows = FALSE,
        angle_col = 90,
        main = pathway)
}

pathways <- grep(colnames(sc_gsea_metadata), pattern = metadata_to_plot, value = TRUE)
for (pathway in pathways) {
    make_heatmap_per_pathway(sc_gsea_metadata = sc_gsea_metadata,
        heatmap_column = heatmap_column,
        heatmap_row = heatmap_row,
        pathway = pathway)
}