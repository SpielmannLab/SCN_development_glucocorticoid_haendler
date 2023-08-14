#!/usr/bin/env Rscript
# To generate a plot of the expression of genes
# conda activate AstizSCN
# Usage:
# figure1_expression_featureplot.R file_sc_obj gene color_low color_high outdir

library(Seurat)
library(ggplot2)
library(dplyr)

# collect arguments
args <- commandArgs(trailingOnly = TRUE)
file_sc_obj <- args[1]
gene <- args[2]
color_low  <- args[3]
color_high  <- args[4]
outdir <- args[5]

assay <- "RNA"
pt_size <- 2
cell_order <- FALSE 
width <- 7
height <- 7

setwd(outdir)

#read suerat object
sc_obj <- readRDS(file_sc_obj)

# Make FeaturePlot
plot <- FeaturePlot(sc_obj,
    features = gene,
    keep.scale = "feature",
    cols = c(color_low, color_high),
    pt.size = pt_size,
    order = cell_order
) &
    theme(legend.position="none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = color_low),
        line = element_blank(),
        plot.title = element_text(family="sans", face="plain", hjust=0.5, size=15),
        aspect.ratio = 1)

filename = paste0(gsub(file_sc_obj,
    pattern = ".rds$",
    replacement = paste0("_FeaturePlot_", gene, "_colors", color_low, "2", color_high, ".tiff"))) %>%
    basename()
ggsave(plot = plot, filename = filename, width = width, height = height, units = "in", dpi = 300)
