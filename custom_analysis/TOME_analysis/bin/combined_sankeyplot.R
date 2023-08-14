#!/usr/bin/env Rscript
#############################
#### Section2_trajectory ####
#############################

#### The scripts used for creating cellular trajectories

##########################
#### Step2_connection_sub step 3 Make Sankey Plots ####
##########################

library(htmlwidgets)
# library(webshot)
library(networkD3)
library(dplyr)
library(ggplot2)

# Get command arguments to get the k value for knn
args <- commandArgs(trailingOnly=TRUE)
k <- as.numeric(args[1])
min_link_strength <- as.numeric(args[2])

# ***** Combined sankey from UMAP ****
files <- list.files(pattern=paste0("k", k, "_umap_sankey_links.rds"))
sankey_links <- lapply(files, readRDS)
combined_sankey_links <- do.call(args=sankey_links, what=rbind)

# Get a distribution of like_strengths and then apply filter
plot <- ggplot(combined_sankey_links, aes(x = value)) + 
    geom_density(fill = "grey", color = "grey") +
    theme_classic()
filename <- paste0("sankey_link_strength_distribution_from_umap_",k,".pdf")
ggsave(plot = plot, filename = filename, width = 7, height = 5)
combined_sankey_links <- filter(combined_sankey_links, value >= min_link_strength)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
name=c(as.character(combined_sankey_links$source), as.character(combined_sankey_links$target)) %>%
  unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
combined_sankey_links$IDsource <- match(combined_sankey_links$source, nodes$name)-1
combined_sankey_links$IDtarget <- match(combined_sankey_links$target, nodes$name)-1

p <- sankeyNetwork(Links = combined_sankey_links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   sinksRight=FALSE, fontSize=20)

filename <- paste0("sankey_from_umap_",k,".html")
saveWidget(p, file=filename)

# ***** Combined sankey from PCA ****
files <- list.files(pattern=paste0("k",k,"_pca_sankey_links.rds"))
sankey_links <- lapply(files, readRDS)
combined_sankey_links <- do.call(args=sankey_links, what=rbind)

# Get a distribution of like_strengths and then apply filter
plot <- ggplot(combined_sankey_links, aes(x = value)) + 
    geom_density(fill = "grey", color = "grey") +
    theme_classic()
filename <- paste0("sankey_link_strength_distribution_from_umap_",k,".pdf")
ggsave(plot = plot, filename = filename, width = 7, height = 5)
combined_sankey_links <- filter(combined_sankey_links, value >= min_link_strength)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
name=c(as.character(combined_sankey_links$source), as.character(combined_sankey_links$target)) %>%
  unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
combined_sankey_links$IDsource <- match(combined_sankey_links$source, nodes$name)-1
combined_sankey_links$IDtarget <- match(combined_sankey_links$target, nodes$name)-1

p <- sankeyNetwork(Links = combined_sankey_links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   sinksRight=FALSE, fontSize=20)

filename <- paste0("sankey_from_pca_",k,".html")
saveWidget(p, file=filename)

# pca_files <- list.files(pattern=paste0(k,"_umap_sankey_links.rds")

# Conversion of html to pdf
# webshot::webshot(filename, gsub(filename, pattern=".html", replacement=".pdf"))
# filename <- paste0("combined_sankey_from_umap_",k_neigh,".rds")
# saveRDS(combined_links, file=filename)
