#!/usr/bin/env Rscript
#Do it in NICHES conda environment
doc <- "Perform differential signalling analysis between the conditions"
print(doc)
rm(doc)

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# --- Parameters: Read variables from the environment and parse as needed
group_by <- Sys.getenv("group_by") %>%
    paste0(".Sending") #Add this, because NICHES does it
signalling_npcs <- Sys.getenv("signalling_npcs")
min_feat_cell2cell <- Sys.getenv("min_feat_cell2cell")

# Print all the variabes
for (i in ls()) if (i != "i") message("The value of ", i, " is: ", get(i))

# ----- Read in the cell signalling objects and merge them
cell2cell_obj1 <- readRDS("cell2cell_obj1.rds")
cell2cell_obj2 <- readRDS("cell2cell_obj2.rds")

cell2cell_obj <- merge(cell2cell_obj1, cell2cell_obj2)

# Do QC to visualize the amount of gene-gene pairs per cell and remove those cells with very little signalling
vln_plot <- VlnPlot(cell2cell_obj,
                    features = 'nFeature_CellToCell',
                    group.by = group_by,
                    pt.size = 0.1,
                    log = TRUE
)
ggsave(vln_plot, filename = "vln_plot_nFeature_cell2cell.pdf", width = 10, height = 7)

cell2cell_obj <- subset(cell2cell_obj,
                        nFeature_CellToCell > min_feat_cell2cell)

# Do basic Seurat analysis on the merged and cleaned up cell2cell signalling data
cell2cell_obj <- cell2cell_obj %>%
    ScaleData(vars.to.regress = "nFeature_CellToCell") %>%
    RunPCA(features = rownames(.)) %>%
    RunUMAP(dims = 1:signalling_npcs)

# This ensures that the idents are sorted properly
Idents(cell2cell_obj) <- as.factor(cell2cell_obj$VectorType)

# Save the Seurat signalling object
saveRDS(cell2cell_obj, file = "combined_cell_signalling_obj.rds")

# Plotting
# First add metadata to plot just by sending or receiving VectorType
cell2cell_obj@meta.data <- cell2cell_obj@meta.data %>%
    mutate(VectorType_sending = gsub(VectorType, pattern = "—.*", replacement = "")) %>%
    mutate(VectorType_receiving = gsub(VectorType, pattern = ".*—", replacement = ""))

sig_plot <- DimPlot(cell2cell_obj,
                    group.by = c("VectorType","VectorType_sending","VectorType_receiving"), 
                    split.by = group_by,
                    shuffle = TRUE,
                    ncol = 1,
                    label = FALSE) &
    theme(aspect.ratio = 1,
         legend.position = "right")
ggsave(sig_plot, filename = "umap_plot.pdf", width = 25, height = 20)

# Split each pairing to get DE genes for each pairing
markers_list <- list()
for (signalling_pair in levels(cell2cell_obj)) {
    message("Finding DE for: ", signalling_pair)
    cell2cell_obj_sub <- subset(cell2cell_obj, idents = signalling_pair)

    # metadata columns after NICHES have ".Sending/.Receiving/.Joint" suffix.
    Idents(cell2cell_obj_sub) <- cell2cell_obj_sub[[group_by]]
    markers <- FindAllMarkers(cell2cell_obj_sub, only.pos = TRUE)

    # Add a column to denote which cluster the markers correspond to as these tables will be concatenated
    # But only if markers were not found
    if(dim(markers)[1] != 0){
        markers$VectorType <- signalling_pair
        markers_list[[signalling_pair]] <- markers
    }
}

# Combined the DE genes calculated between the two groups for each cell2cell cluster
markers_combined <- do.call(markers_list,
                            what = rbind)

write.table(markers_combined,
            file = paste0("DiffSig_between_", group_by, ".tsv"),
            sep = "\t")

# Make Violin plot
markers <- markers_combined %>%
    mutate(sending = gsub(VectorType, pattern = "—.*", replacement = "")) %>%
    mutate(receiving = gsub(VectorType, pattern = ".*—", replacement = "")) %>%
    mutate(cluster = as.factor(cluster))

# Flip the sign
markers$avg_log2FC[markers$cluster == levels(markers$cluster)[2]] <- -markers$avg_log2FC[markers$cluster == levels(markers$cluster)[2]]

plot_sending <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sending)) +
    geom_point() +
    xlab(paste0("<---", levels(markers$cluster)[2], " | avg_log2FC | ", levels(markers$cluster)[1], "--->"))
plot_receiving <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = receiving)) +
    geom_point() +
    xlab(paste0("<---", levels(markers$cluster)[2], " | avg_log2FC | ", levels(markers$cluster)[1], "--->"))
plot_VectorType <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = VectorType)) +
    geom_point() +
    xlab(paste0("<---", levels(markers$cluster)[2], " | avg_log2FC | ", levels(markers$cluster)[1], "--->"))

plot_sending_receiving <- plot_grid(plot_sending, plot_receiving, ncol = 2)
plot <- plot_grid(plot_sending_receiving, plot_VectorType, ncol = 1)
ggsave(plot,
       filename = paste0("Volcano_of_DiffSig_between_", group_by, ".pdf"),
       width = 15, height = 10)
