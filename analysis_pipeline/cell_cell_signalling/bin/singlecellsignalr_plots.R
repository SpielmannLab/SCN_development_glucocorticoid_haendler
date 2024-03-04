#!/usr/bin/env Rscript
# Use MegaBundle
" Make plots (e.g., heatmap of the LR strength for each LR pair and cell-cell pair)

Usage: make_plots.R --file_signal=<file> --min_LRscore=<val> --sampling=<value>

Options:
  -h --help			Show this screen.
  --file_signal=<file>		Seurat Object to be downsampled
  --min_LRscore=<value>   value between 0 and 1. 0.6 is a good start
  --sampling=<value>    For plotting purposes, number of LR pairs to plot
" -> doc

suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyr))
suppressMessages(library(docopt))

# --- Parameters: Read
arguments <- docopt(doc, quoted_args=TRUE)
file_signal <- arguments$file_signal
min_LRscore <- as.numeric(arguments$min_LRscore)
sampling <- as.numeric(arguments$sampling)

message("file_signal: ", file_signal)
message("min_LRscore: ", min_LRscore)
message("sampling: ", sampling)

## Read the data
signal <- readRDS(file_signal)

# Lets generate a heatmap
signal <- lapply(X = signal, FUN = function(x){
    x$sender <- names(x)[1]
    x$receiver <- names(x)[2]
    names(x) <- c("Ligand", "Receptor", "interaction_type", "LRscore", "sender", "receiver")
    x
})
names(signal) <- ""

signal_for_heatmap <- do.call(what = rbind, args = signal) %>%
    filter(LRscore > min_LRscore) %>%
    mutate(LR = paste0(Ligand, "-", Receptor),
        cell_cell = paste0(sender, " -> ", receiver)) %>%
    select(cell_cell, LR, interaction_type, LRscore)

matrix_for_heatmap <- signal_for_heatmap %>%
    select(cell_cell, LR, LRscore) %>%
    distinct(cell_cell, LR, .keep_all = TRUE) %>%
    pivot_wider(names_from = "cell_cell", values_from = "LRscore", values_fill = 0) %>%
    tibble::column_to_rownames(var = "LR") %>%
    as.matrix()

# sample down for visualization
if (sampling < NROW(matrix_for_heatmap) & sampling > 0) {
    sampled_rows <- sample(x = seq.int(from = 1, to = NROW(matrix_for_heatmap), by = 1), size = sampling)
} else {
   sampled_rows <- seq.int(from = 1, to = NROW(matrix_for_heatmap), by = 1)
}

matrix_for_heatmap <- matrix_for_heatmap[sampled_rows, ] %>%
    t() #transpose, so the long cell names are the rownames

file_heatmap <- gsub(file_signal, pattern=".rds", replacement="_heatmap.pdf")
pheatmap(matrix_for_heatmap,
    cluster_cols = FALSE,
    treeheight_col = 20,
    cluster_rows = TRUE,
    treeheight_row = 20,
    filename = file_heatmap,
    cellwidth = 10,
    cellheight = 10)

quit("no")
signal_difference <- inner_join(x=pnd2, y=pnd30, by=c("cell_cell","LR"),suffix=c(".1",".2")) %>% mutate(delta_LRscore = LRscore.1-LRscore.2) %>% select(cell_cell, LR, delta_LRscore)
matrix_for_heatmap <- signal_difference %>% pivot_wider(names_from = "cell_cell", values_from = "delta_LRscore", values_fill = NA) %>%
    tibble::column_to_rownames(var = "LR") %>%
    as.matrix() %>%
    t()

min_max <-  max(abs(c(min(matrix_for_heatmap, na.rm = TRUE), max(matrix_for_heatmap, na.rm = TRUE))))
pheatmap(matrix_for_heatmap,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    breaks = seq(from = -min_max,
        to = +min_max,
        length.out = 101),
    cluster_cols = FALSE, 
    cluster_rows = FALSE, 
    main = "LR_pnd2 - LR_pnd30", 
    cellwidth=40, cellheight=40, 
    display_numbers = TRUE, 
    filename = "pnd2-pnd30_LR.pdf", 
    fontsize = 10)