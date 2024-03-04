#!/usr/bin/env Rscript
#Do it in MegaBundle/scVelocity conda environment
doc <- "Write Seurat Single cell Object as h5ad scanpy file for Liana

Usage: rds_h5ad.R --scobject=<file> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --scobject=<file>    *.rds file LoomObject.

"

# --- Dependencies
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratDisk))
suppressMessages(library(docopt))

arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

# ---------------
# --- Parameters
# ---------------
sc_file <- arguments$scobject
samplename <- gsub(pattern = ".rds", replacement = "", sc_file)
sc <- readRDS(sc_file)
print(sc)

SaveH5Seurat(sc, filename = paste0(samplename, ".h5Seurat"))
Convert(paste0(samplename, ".h5Seurat"), dest = "h5ad")
message(paste0("The h5ad object ", samplename, ".h5ad ", "has been created."))
