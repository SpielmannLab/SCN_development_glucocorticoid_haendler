#!/usr/bin/env Rscript
#Do it in NICHES conda environment
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

# --- Parameters: Read
arguments <- docopt(doc, quoted_args=TRUE)
assay <- arguments$assay

message("assay: ", assay)

# ----- Set default assay
sc_obj <- readRDS("sc_obj.rds")
DefaultAssay(sc_obj) <- assay

# ----- Print to screen the default assay as a check
Assays(sc_obj)

# ----- Save the seurat object with default assay
saveRDS(sc_obj, file = "sc_obj_default_assay_set.rds")
