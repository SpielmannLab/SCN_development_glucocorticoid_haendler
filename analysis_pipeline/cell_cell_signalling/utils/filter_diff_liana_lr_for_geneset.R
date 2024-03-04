#!/usr/bin/env Rscript
# ----- Load libraries ------#
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
# Define parameters

# get all the arguments and files
file_diff_liana_scores <- "/Users/sreenivasan/Documents/Works/AstizSCN/cell_cell_signalling_analysis/cell_signalling_output/morris_equal_cells/diff_liana/morris_scn_scores/Dysregulated_LRs_unfiltered.tsv" 
tileplot_width <- 12
tileplot_height <- 12
celltypes_to_plot = c("01 SCN neurons", "03 Astrocytes", "04 Ependymal cells")
file_collecTRI_regulon <- "/Users/sreenivasan/Documents/Works/AstizSCN/cell_cell_signalling_analysis/cell_signalling_output/liana_custom_lr/collectTRI_Nr3c1_regulons.txt" 

# Read in the files
lr_res_alt = read.table(file_diff_liana_scores, sep = "\t", header = TRUE, row.names = 1)
collecTRI_regulon <- read.table(file_collecTRI_regulon, header = TRUE, row.names = 1) |>
  select(-source)

# Define a file to keep logs
message("Number of genes activated/repressed (1/-1) by Nr3c1:")
table(collecTRI_regulon$mor)

# Filter ligands or receptors that are significantly regulated
# in the correct direction
lr_res_alt_filtered <- lr_res_alt |>
    dplyr::left_join(collecTRI_regulon, by = c("ligand" = "target"), suffix <- c("", ".ligand"), keep = FALSE) |>
    dplyr::left_join(collecTRI_regulon, by = c("receptor" = "target"), suffix <- c("", ".ligand"), keep = FALSE, suffix = c(".ligand", ".receptor")) |>
    filter(!is.na(mor.ligand) | !is.na(mor.receptor)) |>
    dplyr::filter(((receptor_stat*mor.receptor > 0) & (receptor_padj < 0.05)) | ((ligand_stat*mor.ligand > 0) & (ligand_padj < 0.05))) |>
    select(-mor.ligand, -mor.receptor)

message("Number of LR pairs in which either the ligand or the receptor is significantly differently modulated by Nr3c1 in the right direction):")
NROW(lr_res_alt_filtered)

# Filter based on source and targets cell types requested to plot
lr_res_alt_filtered <- lr_res_alt_filtered |>
      dplyr::filter((source %in% celltypes_to_plot) & (target %in% celltypes_to_plot))
message("Number of LR pairs remaining in the cell types to plot")
NROW(lr_res_alt_filtered)

# Save the filtered lr pairs for plotting in python
write.table(lr_res_alt_filtered,
    file = "/Users/sreenivasan/Documents/Works/AstizSCN/cell_cell_signalling_analysis/cell_signalling_output/liana_custom_lr/sig_lr_in_collecTRI_Nr3c1_regulon.tsv",
    sep = "\t",
    row.names = FALSE)
