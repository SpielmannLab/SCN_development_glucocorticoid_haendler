#!/usr/bin/env Rscript
# Do it in NICHES conda environment
doc <- " Perform the basic NICHES analysis to convert single-cell object to single-cell signalling object. Will use the default assay

Usage: do_niches.R

It also requires the following variables as environmental variables:
    species (human|mouse)
    LRdatabase (omnipath|phantom)
    cell_types_column (metadata column name containing cell type annotation)
    min_cells_p_gene (number)
    min_cells_p_cluster (number)
    mode_of_analysis (CellToCell|CellToSystem|SystemToCell)
    signalling_npcs (Number of pcs for dimreduction of cellsignalling object)
"
print(doc)
rm(doc)
# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(NICHES))
suppressMessages(library(ggplot2))

# --- Parameters: Read variables from the environment and parse as needed
species <- Sys.getenv("species")
lr_database <- Sys.getenv("LRdatabase")
cell_types_column <- Sys.getenv("cell_types_column")
min_cells_p_gene <- Sys.getenv("min_cells_p_gene") %>%
  as.integer()
min_cells_p_cluster <- Sys.getenv("min_cells_p_cluster") %>%
  as.integer()
mode_of_analysis <- Sys.getenv("mode_of_analysis")
mode_of_analysis <- mode_of_analysis %>%
  gsub(pattern = "(\\[|\\])", replacement = "") %>%
  strsplit(",") %>%
  unlist()
signalling_npcs <- Sys.getenv("signalling_npcs")

# Print all the variables
for (i in ls()) if (i != "i") message("The value of ", i, " is: ", get(i))

# ----- Check which assay to use using the DefaultAssay
sc_obj <- readRDS("sc_obj.rds")
assay <- DefaultAssay(sc_obj)
message("The assay, to be used based on DefaultAssay: ", assay)

# ----- Perform NICHES
# This prepping step needs to be done separately, because the slot Count in the alra assay is undefined
Idents(sc_obj) <- sc_obj[[cell_types_column]]
sc_obj <- prepSeurat(
  object = sc_obj,
  assay = "RNA",
  min.cells.per.ident = min_cells_p_cluster,
  min.cells.per.gene = min_cells_p_gene
)
# Set the default assay back, because prepSeurat sets the default assay to RNA
DefaultAssay(sc_obj) <- assay

sc_sig_obj <- RunNICHES(
  object = sc_obj,
  assay = assay,
  LR.database = lr_database,
  species = species,
  cell_types = cell_types_column,
  CellToCell = ("CellToCell" %in% mode_of_analysis),
  CellToSystem = ("CellToSystem" %in% mode_of_analysis),
  SystemToCell = ("SystemToCell" %in% mode_of_analysis)
)

# ----- Do further plooting of CellToCell signalling if asked
if (!("CellToCell" %in% mode_of_analysis)) {
  quit("no")
}

cell2cell_obj <- sc_sig_obj$CellToCell %>%
  ScaleData() %>%
  RunPCA(features = rownames(.)) %>%
  RunUMAP(dims = 1:signalling_npcs)

# This ensures that the idents are sorted properly
Idents(cell2cell_obj) <- as.factor(cell2cell_obj$VectorType)

# ----- Save the seurat object containing the cell2cell signalling data
saveRDS(cell2cell_obj, file = "CellToCell_obj.rds")

eb_plot <- ElbowPlot(cell2cell_obj, ndims = 50)
ggsave(eb_plot, filename = "elbow_plot.pdf")

sig_plot <- DimPlot(cell2cell_obj, group.by = "VectorType", label = FALSE) +
  theme(aspect.ratio = 1)
ggsave(sig_plot, filename = "umap_plot.pdf", width = 10, height = 10)

# Find markers
markers <- FindAllMarkers(cell2cell_obj, min.pct = 0.25, only.pos = T, test.use = "roc")
GOI_niche <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, myAUC)

marker_plot <- DoHeatmap(cell2cell_obj, features = unique(GOI_niche$gene)) +
  scale_fill_gradientn(colors = c("grey", "white", "blue"))
ggsave(marker_plot, filename = "cluster_wise_marker_plot.pdf", width = 10, height = 10)
