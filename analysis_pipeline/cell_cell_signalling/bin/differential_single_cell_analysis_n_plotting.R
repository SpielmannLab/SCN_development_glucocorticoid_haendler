#!/usr/bin/env Rscript
# conda activate NICHES

doc <- "Differential analysis of the input Seurat object"
print(doc)
rm(doc)

# --- load libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))

# --- Parameters: Read variables from the environment and parse as needed
group_by <- Sys.getenv("group_by")
cell_types_column <- Sys.getenv("cell_types_column")

# ----- Read in the single-cell seurat object
sc_obj <- readRDS("sc_obj.rds")

# Make Plot of sequencing depth
plot <- VlnPlot(sc_obj, group.by = cell_types_column, split.by = group_by, features = c("nCount_RNA", "nFeature_RNA")) &
    theme(legend.position = "bottom")
ggsave(plot, filename = "VlnPlot_nCount_nFeature.pdf", width = 7, height = 4)

# Differential gene expression analysis per cell type
# Split each cell type to get DE genes for each cell type
Idents(sc_obj) <- sc_obj[[cell_types_column]]
markers_list <- list()
for (cell_type in levels(sc_obj)) {
    message("Finding DE for: ", cell_type)
    sc_obj_sub <- subset(sc_obj, idents = cell_type)

    # metadata columns after NICHES have ".Sending/.Receiving/.Joint" suffix.
    Idents(sc_obj_sub) <- sc_obj_sub[[group_by]]
    markers <- FindAllMarkers(sc_obj_sub, only.pos = TRUE)

    # Add a column to denote which cluster the markers correspond to as these tables will be concatenated
    # But only if markers were not found
    if(dim(markers)[1] != 0){
        markers$cell_type <- cell_type
        markers_list[[cell_type]] <- markers
    }
}

# Combined the DE genes calculated between the two groups for each cell type
markers_combined <- do.call(markers_list,
                            what = rbind)

write.table(markers_combined,
            file = paste0("DE_between", group_by, ".tsv"),
            sep = "\t")

# Make Violin plots
markers <- markers_combined %>%
    mutate(cluster = as.factor(cluster))
# Flip the sign
markers$avg_log2FC[markers$cluster == levels(markers$cluster)[2]] <- -markers$avg_log2FC[markers$cluster == levels(markers$cluster)[2]]

plot <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = cell_type)) + geom_point()+
    xlab(paste0("<---", levels(markers$cluster)[2], " | avg_log2FC | ", levels(markers$cluster)[1], "--->"))
ggsave(plot, filename = "volcano_of_diff_expression.pdf", width = 8, height = 6)
