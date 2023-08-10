#!/usr/bin/env Rscript
# Call from conda-environment MegaBundle
" Make heatmap plots based on the GSEA analysis results

Usage: gsea_heatmap.R --file_sc_gsea_metadata=<file> --metadata_to_plot=<value> --groupby=<value>

Options:
  -h --help			Show this screen.
  --file_sc_gsea_metadata=<file>		RDS file containing sc_obj metadata and ssGSEA scores
  --metadata_to_plot=<value>  Use <ESCAPE_>. This has been appendend to all ssGSEA enrichment score colnames
  --groupby=<value>   Which metadata column of Seurat object to group for statistical testing and violin plotting
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
message("groupby: ", groupby)

# --- Read files
message("Loading Data")
sc_gsea_metadata <- readRDS(file_sc_gsea_metadata)

# get median gsea scores segregated by annotations for each of the Gene sets
medians_all <- sc_gsea_metadata %>%
    group_by(across(all_of(groupby))) %>%
    dplyr::summarise(across(starts_with(metadata_to_plot), mean)) %>%
    tibble::column_to_rownames(var = groupby)

#define a function to drop names
drop_escape_in_names <- function(x){
    names <- colnames(x)
    colnames(x) <- gsub(names,pattern="ESCAPE_", replacement="") %>%
        gsub(,pattern="_", replacement=" ")
    return(x)
}

medians <- medians_all %>%
    as.matrix() %>% # convert to matrix for pheatmap
    drop_escape_in_names() # cleanup names - see names above

filename = paste0("Heatmap_",gsub(file_sc_gsea_metadata, pattern=".rds", replacement=paste0("_", groupby, ".pdf")))
pheatmap(medians, cellheight = 20, cellwidth = 20, filename = filename, cluster_rows = FALSE, angle_col = 90)