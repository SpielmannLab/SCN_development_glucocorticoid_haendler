#!/usr/bin/env Rscript

" Run the SingleCellSignalR pipeline

# Follows vignette in the BioConductor page as well as github page
# https://www.bioconductor.org/packages/release/bioc/vignettes/SingleCellSignalR/inst/doc/UsersGuide.html
# https://github.com/SCA-IRCM/Demo

Usage: singlecellsignalr.R --file_matrix=<file> --file_clusters=<file> --file_annotations=<file> --species=<value> --min_logFC=<value>

Options:
  -h --help			Show this screen.
  --file_matrix=<file>		Input file exported using seurat_to_matrix
  --file_clusters=<file>		Input file exported using seurat_to_matrix
  --file_annotations=<file>		Input file exported using seurat_to_matrix
  --species=<value> 'mus musculus' or 'homo sapiens'
  --min_logFC=<value>   The log fold-change threshold for differentially expressed genes. Use 2 as a start

" -> doc

# --- load libraries
suppressMessages(library(dplyr))
suppressMessages(library(SingleCellSignalR))
suppressMessages(library(docopt))

# --- Parameters: Read
arguments <- docopt(doc, quoted_args=TRUE)
file_matrix <- arguments$file_matrix
file_clusters <- arguments$file_clusters
file_annotations <- arguments$file_annotations
species <- arguments$species
min_logFC <- arguments$min_logFC

message("file_matrix: ", file_matrix)
message("file_clusters: ", file_clusters)
message("file_annotations: ", file_annotations)
message("species: ", species)
message("min_logFC: ", min_logFC)

## Read the data
matrix <- readRDS(file_matrix)
clusters <- readRDS(file_clusters)
annotations <- readRDS(file_annotations)

annotations = annotations %>%
    gsub(annotations, pattern="\\W", replacement = "_")

all.genes <- rownames(matrix)

# Analyse the interactions
signal <- cell_signaling(data = matrix, 
    genes = all.genes,
    cluster = clusters,
    c.names = annotations,
    species = species,
    logFC = logFC,
    write = FALSE)

# Save the interaction values to yaml
file_signal <- gsub(file_matrix, pattern=".rds", replacement="_LRinteractions.rds")
saveRDS(signal, file = file_signal)

# Visualization and plot export
# A summary of the interactions between cell clusters can be output in the form of a chord diagram
plot_title <- basename(file_matrix)
visualize_interactions(signal, write.out = TRUE, method = plot_title)