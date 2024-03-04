#!/usr/bin/env Rscript
" Export count or data matrix, clusternumber, and annotation from Seurat object

Usage: seurat_to_matrix.R --file_sc_obj=<file> --slot=<value> --assay=<value> --cluster_column=<value> --annotation_column=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		Seurat Object to be downsampled
  --assay=<value>   'RNA', 'SCT', 'integrated'
  --slot=<value>  'data' for normalized, 'counts' for raw counts.
  --cluster_column=<value>  Column name of the cluster
  --annotation_column=<value>   Column name containing annotations for the clusters
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))

# --- Parameters: Read
arguments <- docopt(doc, quoted_args=TRUE)
file_sc_obj <- arguments$file_sc_obj
assay <- arguments$assay
slot <- arguments$slot
cluster_column <- arguments$cluster_column
annotation_column <- arguments$annotation_column

message("file_sc_obj: ", file_sc_obj)
message("assay: ", assay)
message("slot: ", slot)
message("cluster_column: ", cluster_column)
message("annotation_column: ", annotation_column)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

if (slot == "data") {
    matrix <- data.frame(sc_obj[[assay]]@data)
} else if (slot == "count") {
    matrix <- data.frame(sc_obj[[assay]]@counts)
}

# Save count matrix
filename=gsub(file_sc_obj, pattern=".rds", replacement="_matrix.rds")
saveRDS(matrix, file = filename)

# save clusters and cluster annotation as a single data.frame
clusters <- as.numeric(sc_obj@meta.data[, cluster_column])
# Since sometimes in subsetted data, not all clusters are present, rename them from 0 onwards:
clusters <- clusters %>%
    as.factor() %>%
    plyr::mapvalues(from=levels(.), seq.int(from = 0, length.out = length(levels(.)))) %>%
    as.numeric(levels(.))[.]
file_clusters <- gsub(file_sc_obj, pattern=".rds", replacement="_clusters.rds")
saveRDS(clusters, file_clusters)

if (!is.null(annotation_column) | !(annotation_column %in% names(sc_obj@meta.data))) {
    annotations <- sc_obj@meta.data %>%
        select(all_of(c(cluster_column, annotation_column))) %>%
        distinct(across(cluster_column),.keep_all = TRUE) %>%
        arrange(across(all_of(cluster_column))) %>%
        pull() %>%
        as.character()
    file_annotations <- gsub(file_sc_obj, pattern=".rds", replacement="_annotations.rds")
    saveRDS(annotations, file_annotations)
}
