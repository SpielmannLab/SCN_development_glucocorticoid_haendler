#############################
#### Section2_trajectory ####
#############################

#### The scripts used for creating cellular trajectories

##########################
#### Step2_2_connection ####
##########################

#### running Knn to find ancestor

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(htmlwidgets)
library(plotly)
library(monocle3)
library(FNN)
source("functions.R")

# Get command arguments to get the k value for knn
args <- commandArgs(trailingOnly=TRUE)
k <- as.numeric(args[1])

# extract just the required metadata
sc_obj = readRDS("file_sc_obj.rds")
anno <- sc_obj@meta.data[,c("orig.ident", "Anno", "day")]
anno$stage = anno$orig.ident

# Get the umap embeddings
emb <- readRDS("file_umap_coords.rds")
if(sum(!(rownames(anno) %in% rownames(emb)))){
    print("Error!")
    print(xxx)
}
anno = anno[rownames(emb),] # This is not necessary

res = createLineage_Knn(emb, anno, k_neigh = k, permute=TRUE) #### createLineage_Knn with permutation
filename <- paste(unique(anno$orig.ident), collapse="_") %>%
  paste0("_Knn_k",k,"_umap_lineage.rds") %>%
  tolower()
saveRDS(res, file=filename)

# save proportion edge weight over 0.2
edge_weights <- as.vector(as.matrix(res))
edge_weight_prop <- sum(edge_weights>0.2)/length(edge_weights)
filename <- paste(unique(anno$orig.ident), collapse="_") %>%
  paste0("_Knn_k",k,"_edge_weight_prop.txt") %>%
  tolower()
write.table(x=edge_weight_prop, file=filename, row.names=FALSE, col.names="Fraction Edge Weight > 2")

# Create links for sankey plot
sankey_links <- make_sankey_links(knn_matrix=res)
filename <- paste(unique(anno$orig.ident), collapse="_") %>%
  paste0("_Knn_k",k,"_umap_sankey_links.rds") %>%
  tolower()
saveRDS(sankey_links, file=filename)

#### repeat this approach but using 30 PCs instead to testing if results are robust
# Get the pca embeddings
emb = readRDS("file_pca_coords.rds")
if(sum(!(rownames(anno) %in% rownames(emb)))){
    print("Error!")
    print(xxx)
}
anno = anno[rownames(emb),]

res = createLineage_Knn(emb, anno, replication_times = 100, k_neigh= k) #### createLineage_Knn function was in help_code.R
sankey_links <- make_sankey_links(knn_matrix=res)

filename <- paste(unique(anno$orig.ident), collapse="_") %>%
  paste0("_Knn_k",k,"_pca_lineage.rds") %>%
  tolower()
saveRDS(res, file=filename)
filename <- paste(unique(anno$orig.ident), collapse="_") %>%
  paste0("_Knn_k",k,"_pca_sankey_links.rds") %>%
  tolower()
saveRDS(sankey_links, file=filename)
