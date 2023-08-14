#!/usr/bin/env Rscript
# Call from conda-environment MegaBundle
" Downsample a seurat object based on the desired cell fraction. Output is an RDS file containing the downsampled Seurat object.

Usage: rscript.R --file_sc_obj=<file> --downsample=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		Seurat Object to be downsampled
  --downsample=<value>  If the value is <=1, treated as a fraction to downsample by. If >1, then dowsamples to this number of cells
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))

# --- Parameters: Read
arguments <- docopt(doc, quoted_args=TRUE)
file_sc_obj <- arguments$file_sc_obj
downsample <- as.numeric(arguments$downsample)

message("file_sc_obj: ", file_sc_obj)
message("downsample: ", downsample)

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

# Downsample using ths subset function in SeuratObject
if (downsample <= 1){
  no_cells_to_subset <- length(Cells(sc_obj))*downsample
} else {
  no_cells_to_subset = downsample
}

set.seed(111)
cells_to_subset <- sample(Cells(sc_obj), size=no_cells_to_subset)
sc_obj_downsampled <- subset(sc_obj, cells=cells_to_subset)

filename=gsub(file_sc_obj, pattern=".rds", replacement="_downsampled.rds")
saveRDS(sc_obj_downsampled, file=filename)
