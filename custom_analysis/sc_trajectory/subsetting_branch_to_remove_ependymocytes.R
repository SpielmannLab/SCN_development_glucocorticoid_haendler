# This script was used on a local computer with GUI. It was used to remove the ependymocytes branch (bottom left) manually. It was done using the function "choose_graph_segments".

suppressMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

setwd("/Users/sreenivasan/Documents/Publications/2023 Glucocorticoid/output/sc_trajectory/astrocytes_by_annotation/local_processing")

file_cds <- "/Users/sreenivasan/Documents/Publications/2023 Glucocorticoid/output/sc_trajectory/astrocytes_by_annotation/trajectory/annotated_annotations_subsetted_noextraSCN_annotations_subsetted_astrocytes_cds_traj.rds"
file_sc_obj <- "/Users/sreenivasan/Documents/Publications/2023 Glucocorticoid/output/sc_trajectory/astrocytes_by_annotation/annotated_annotations_subsetted_noextraSCN_annotations_subsetted_astrocytes.rds"
assay <- "RNA"
genes <- c("Nr3c1", "Nr3c2")
gex_pt_size <- 2
gex_width <- 7
gex_height <- 7
filename_prefix <- "astroytes_by_annotation_without_ependymocytes"

# Open the parent seurat object
sc_obj <- readRDS(file_sc_obj)
DefaultAssay(sc_obj) <- assay

cds <- readRDS(file_cds)
cds_sub <- choose_graph_segments(cds) # Choose starting node based on the automated start node in the cds obj. Then choose all the 4 astrocytes branch points as end nodes. 

cds <- cds_sub


# Redo UMAP, clustering, and trajectory construction on the subsetted dataset
npcs <- 5
clustering_k <- 10
seperate_trajectory_by_partition <- FALSE # logical. Whether trajectories should be disconnected between partitions
close_loop <- FALSE # logical. Whether circular trajectories are allowed
group_bys <- "developmental_timepoint"
root_metadata_key <- "developmental_timepoint"
root_metadata_val <- "GD17.5"
keep_embeddings <- FALSE
root_node <- NULL
pt_size <- 2
width <- 12
height <- 6

cds <- preprocess_cds(cds, num_dim = npcs) #does normalization, pca etc.

# Whether to create a new UMAP embedding or to use the one from Seurat
if (keep_embeddings) {
    # replace the pca, as the user wants
    reducedDim(cds, type = "PCA") <- Embeddings(sc_obj, reduction = "pca")
    reducedDim(cds, type = "UMAP") <- Embeddings(sc_obj, reduction = "umap")
} else {
    method_to_use <- "PCA" #alternative is "Aligned"
    seed.use <- 111
    cds <- reduce_dimension(cds, preprocess_method = method_to_use)
}


# Cluster and partition the cells
cds <- cluster_cells(cds, reduction_method = "UMAP", k = clustering_k, random_seed=111)

# Make trajectory
seed.use <- 111
cds <- learn_graph(cds,
    learn_graph_control = list(minimal_branch_len = 20),
    use_partition = seperate_trajectory_by_partition,
    close_loop = close_loop)


get_earliest_principal_node <- function(cds, root_metadata_key, root_metadata_val){
    cell_ids <- which(colData(cds)[, root_metadata_key] == root_metadata_val)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    return(root_pr_nodes)
}
# Fix root_node
if ((root_node %in% c("null", "NULL")) || is.null(root_node) || is.na(root_node)) {
    key_not_avail <- (root_metadata_key %in% c("null", "NULL")) || is.null(root_metadata_key) || is.na(root_metadata_key)
    val_not_avail <- (root_metadata_val %in% c("null", "NULL")) || is.null(root_metadata_val) || is.na(root_metadata_val)
    if (any(key_not_avail, val_not_avail)){
        stop("Neither root_node or root_metadata_key/val pair provided")
    } else {
        root_node <- get_earliest_principal_node(cds, root_metadata_key, root_metadata_val)
    }
}

message("root_node to be used: ", root_node)

# order cells by peudotime
cds <- order_cells(cds, root_pr_nodes = root_node)

# save the cds object
filename <- paste0(filename_prefix, "_pseudotime.rds")
saveRDS(object = cds, file = filename)

# Save monocle plots
plots1 <- lapply(X = group_bys, FUN = function(group_by) {
    plot_cells(cds,
        color_cells_by = group_by,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        graph_label_size = 4,
        cell_size = pt_size) +
        ggtitle(group_by) +
        theme_void() +
        theme(aspect.ratio = 1)
})

plots2 <- lapply(X = c("cluster", "pseudotime"), FUN = function(group_by) {
    plot_cells(cds,
        color_cells_by = group_by,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        graph_label_size = 4,
        cell_size = pt_size) +
        ggtitle(group_by) +
        theme_void() +
        theme(aspect.ratio = 1)
})

#rasterize the plots before saving it to pdf
# suppressPackageStartupMessages(library(ggrastr))
# p <- lapply(p, function(plot){
#   plot <- rasterize(plot, layers="Point", dpi=300)
#   return(plot)
# })

plots <- c(plots1, plots2) # combine all the plots

filename <- paste0(filename_prefix, "_pseudotime",  ".png")
ncol <- length(plots) %>%
    sqrt() %>%
    ceiling()
ggsave(plot <- plot_grid(plotlist = plots, ncol = ncol),
    filename = filename,
    width = width,
    height = height)

# Extract the pseudotime information from cds and add it to Seurat object as new Metadata
sc_obj_subset <- subset(sc_obj, cells = Cells(cds))
if(!identical(Cells(sc_obj_subset), Cells(cds))) stop("The cells in cds obj does not match with that in the seurat object for pseudotime transfer")

sc_obj_subset$pseudotime <- pseudotime(cds)

filename <- paste0(filename_prefix, "_traj_gex_in_pseudotime.rds")
saveRDS(sc_obj_subset, file = filename)

##### Now make plöts based on this Seurat object
gex_pt_size <- 2
gex_width <- 7
gex_height <- 7

group_by <- "developmental_timepoint"
colors_timepoints <- c("zero-expression", "GD17.5", "PND02", "PND10", "PND30")
colors <- c("lightgrey", "#a1dab4","#41b6c4","#2c7fb8","#253494")


for (gene in genes) {
    df_to_plot <- data.frame(expression = as.vector(FetchData(sc_obj_subset, vars = gene, slot = "data")),
            pseudotime = sc_obj_subset$pseudotime,
            group_by = sc_obj_subset$developmental_timepoint)
    colnames(df_to_plot)[1] <- "expression"

    df_to_plot$color_by <- df_to_plot$group_by
    df_to_plot$color_by[df_to_plot$expression == 0] <- "zero-expression"

    plot <- ggplot(df_to_plot, aes(x = pseudotime, y = expression, color = color_by)) +
            geom_point(size = gex_pt_size,
            shape = 16) +
        scale_color_manual(values = colors,
            breaks = colors_timepoints) +
        geom_smooth(formula = y ~ splines::ns(x, df=2),
            se = FALSE,
            color = "black") +
        ylim(0, NA) +
        theme_classic() +
        theme(aspect.ratio = 0.5) 

    filename <- paste0(filename_prefix, "_traj_gex_in_pseudotime_", gene, ".pdf")
    ggsave(plot = plot, filename = filename, width = gex_width, height = gex_height)
}

