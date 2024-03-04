#!/usr/bin/env Rscript
# Do it in NICHES conda environment
"Impute dropout events using ALRA (https://github.com/KlugerLab/ALRA), using SeuratWrappers. Require a normalized matrix, but will check for normalization and perform it if needed.

Usage: impute_by_alra --assay=<value>

Options:
-h --help			Show this screen.
--assay=<value>   Which assay to use for imputing and cell-signalling calculation. Options are <RNA>, <SCT>, <integrated>
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(NICHES))
suppressMessages(library(ggplot2))

# --- Parameters: Read
arguments <- docopt(doc, quoted_args = TRUE)
assay <- arguments$assay

message("assay: ", assay)

# ----- Check if normalized, if not do it
sc_obj <- readRDS("sc_obj.rds")
not_normalized <- identical(
  GetAssayData(sc_obj, assay = assay, slot = "data"),
  GetAssayData(sc_obj, assay = assay, slot = "count")
)
if (not_normalized) {
  warning("The seurat object has not been normalized. Doing lognorm now")
  sc_obj <- NormalizeData(sc_obj)
}

# ----- Perform ALRA
sc_obj <- SeuratWrappers::RunALRA(sc_obj)

# ----- Save imputed seurat object
saveRDS(sc_obj, file = "sc_obj_imputed.rds")
