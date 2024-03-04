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

# Remove scale.data from the object, so that the h5ad object gets the raw counts
sc <- DietSeurat(
  sc,
  counts = TRUE, # so, raw counts save to adata.raw.X
  data = TRUE, # so, log1p counts save to adata.X
  scale.data = FALSE, # set to false, or else will save to adata.X
  features = rownames(sc), # export all genes, not just top highly variable genes
  assays = Assays(sc),
  dimreducs = Reductions(sc),
  graphs = Graphs(sc),
  misc = TRUE
)

print(sc)

SaveH5Seurat(sc, filename = paste0(samplename, ".h5Seurat"))
Convert(paste0(samplename, ".h5Seurat"), dest = "h5ad")
message(paste0("The h5ad object ", samplename, ".h5ad ", "has been created."))
