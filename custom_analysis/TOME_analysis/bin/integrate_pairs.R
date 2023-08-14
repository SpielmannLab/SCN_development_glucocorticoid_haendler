#!/usr/bin/env Rscript
#### The scripts used for creating cellular trajectories

#### Reference: https://satijalab.org/seurat/archive/v3.0/integration.html
#### under "Standard Workflow"

library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2, future.seed=TRUE)
source("functions.R")

### adding both time and batch (if necessary) information
sc_obj_1 = readRDS("file_sc_obj_1.rds") #This was the name given for the files in nextflow
sc_obj_1$day <- "pre"

sc_obj_2 = readRDS("file_sc_obj_2.rds")
sc_obj_2$day <- "nex"

# Add prefix to cell barcodes, so they are unique when merged.
if (length(unique(sc_obj_1$orig.ident))==1 && length(unique(sc_obj_2$orig.ident))==1){
  cell.ids <- c(sc_obj_1$orig.ident[1], sc_obj_2$orig.ident[1]) %>%
    as.character()
}

sc_obj <- merge(sc_obj_1, sc_obj_2, add.cell.ids=cell.ids)

# do standard Seurat integration, since the datasets are small (10x)
# Lot of parameters are defined in the function. Check.
sc_obj.integrated <- doSeuratIntegration(sc_obj)

filename <- paste(unique(sc_obj$orig.ident), collapse="_") %>%
  paste0("_integrated.rds") %>%
  tolower()
saveRDS(sc_obj.integrated, file=filename)

# Save embeddings figure and co-ordinates
filename <- paste(unique(sc_obj.integrated$orig.ident), collapse="_") %>%
  paste0("_tsne.png") %>%
  tolower()
png(filename, 15, 5, units="in", res=300)
print(DimPlot(sc_obj.integrated, group.by=c("orig.ident", "Anno", "ident"), reduction = "tsne", label = TRUE, pt.size = 1, shuffle=TRUE) &
  NoLegend() &
  theme(aspect.ratio=1))
dev.off()

filename <- paste(unique(sc_obj.integrated$orig.ident), collapse="_") %>%
  paste0("_umap.png") %>%
  tolower()
png(filename, 15, 5, units="in", res=300)
print(DimPlot(sc_obj.integrated, group.by=c("orig.ident", "Anno", "ident"), reduction = "umap", label = TRUE, pt.size = 1, shuffle=TRUE) &
  NoLegend() &
  theme(aspect.ratio=1))
dev.off()

emb = data.frame(Embeddings(object = sc_obj.integrated, reduction = "umap"))
filename <- paste(unique(sc_obj$orig.ident), collapse="_") %>%
  paste0("_umap_coords.rds") %>%
  tolower()
saveRDS(emb, file=filename)

emb = data.frame(Embeddings(object = sc_obj.integrated, reduction = "pca"))
filename <- paste(unique(sc_obj$orig.ident), collapse="_") %>%
  paste0("_pca_coords.rds") %>%
  tolower()
saveRDS(emb, file=filename)
