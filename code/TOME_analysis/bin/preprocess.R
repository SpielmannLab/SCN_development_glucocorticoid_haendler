#!/usr/bin/env Rscript
# Run as Rscript preprocess.R annotations_column sex
# This script is used to do do initial bits of preprocessing before calculating linkages.
# Do clustering, Name clusters based on annotations, filter cells based on sex (Xist), find markers, Make lots of plot,

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

#set resolution
args <- commandArgs(trailingOnly = TRUE)
annotations_column <- args[1]
sex <- args[2]

#read file
sc_obj <- readRDS("file_sc_obj.rds")

# Store the requested annotations as Idents
Idents(sc_obj) <- sc_obj[[annotations_column]]
sc_obj$Anno <- paste(sc_obj$orig.ident, Idents(sc_obj), sep="_")

if(sex == "both"){
  message("Keeping cells from both sexes")
} else if(sex == "male"){
  message("Keeping only male cells")
  sc_obj <- subset(sc_obj, subset = sex_Lognorm == 0)
} else if(sex =="female"){
  message("Keeping only female cells")
  sc_obj <- subset(sc_obj, subset = sex_Lognorm > 0)
}

# Save the file
filename <- unique(sc_obj$orig.ident) %>%
  paste0("_processed.rds") %>%
  tolower()
saveRDS(sc_obj, file = filename)

# Save all the generated plots as a single pdf
filename <- unique(sc_obj$orig.ident) %>%
  paste0("_dimplot_phylo_heatmap.pdf") %>%
  tolower()
pdf(filename, width=15, height=1*length(levels(sc_obj)))

sc_obj <- BuildClusterTree(sc_obj,
  features = VariableFeatures(sc_obj),
  reorder = FALSE)
markers <- FindAllMarkers(sc_obj,
  logfc.threshold = 0.25,
  min.pct = 0.2,
  test.use = "wilcox",
  only.pos = TRUE)

# save the markers
filename <- unique(sc_obj$orig.ident) %>%
  paste0("_markers.tsv") %>%
  tolower()
write.table(markers, file = filename, sep = "\t")

top_5_markers <- markers %>%
  group_by(cluster) %>%
  slice_head(n=5)

p1 <- DimPlot(sc_obj, label = TRUE, pt.size = 1) + theme(aspect.ratio = 1)
p2 <- PlotClusterTree(sc_obj)
# Plot only a maximum of 100 cells per group
cells_heatmap <- data.frame(cells=Cells(sc_obj), ident=Idents(sc_obj)) %>%
 group_by(ident) %>%
 slice_sample(n=2) %>%
 .$cells %>%
 head
p3 <- DoHeatmap(
  sc_obj,
  features = top_5_markers$gene,
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  size=3,
  angle=90,
  slot = "scale.data",
  label = TRUE,
  raster = FALSE,
  draw.lines = TRUE,
  combine = TRUE) +
  scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))

# Plot features from Campbell et al and Hajdarovic
features_campbell <- c("Ptgds", "Car2", "Trf", "Plp1", "Mag", "Cldn11","Hbb-bs","Hbb-bt", "Hba-a1", "Rgs5", "Acta2","Igfbp7","Bcas1","Gpr17","Fyn","Pdgfra","Fabp7","Lhfpl3","Bmp4", "Neu4", "Tmem2", "Rras2", "Opalin", "Plekhh1", "Mog","Ctss","C1qb","C1qa","Dcn","Mgp","Apod","Ccdc153","Tmem212","Rarres2","Agt","Kcnd2","Slc1a2","Slc1a3","Scn7a","Cldn10","Col23a1", "Gpr50", "Rax", "Frzb", "Efna5",  "Slit2", "Ch4dl1", "Igfbp2", "Islr", "Fa2h", "Prr5l", "St18", "Col25a1", "6330403K07Rik", "Nnat", "Crym", "Oxt","Pmch","Atp1b1","Vlp","Rgs16","Dlk1","Tac2","Prlr","Tmem35","Gal","Ghrh","Penk","Syt1","Slc18a2","Coch","Npy","Agrp","Crabp1","Ces1d","Epcam","Cyp2f2","Chga","Chgb","Scg2")
p4 <- DotPlot(sc_obj, features=features_campbell) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Gene") +
  ylab("Cluster") +
  ggtitle("Known Markers Campbell 10.1038/nn.4495, Hajdarovic 10.1038/s43587-022-00246-4 and La Manno 10.1038/s41586-018-0414-6") +
  theme(legend.position="bottom", aspect.ratio=length(unique(Idents(sc_obj)))/length(features_campbell) + 0.1)

# Make dotplot of the genes to mark extra-SCN neurons from Hypothalamus, based on Morris et al doi: 10.15252/embj.2021108614
features_extra_scn <-  c("Sst", "Npy", "Sncb", "Pcp4", "Ddit3", "Gal")
p5 <- DotPlot(sc_obj, features=features_extra_scn) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Gene") +
  ylab("Cluster") +
  ggtitle("Expression of extra-SCN (hypothalamic) neuronal markers from Morris et al") +
  theme(legend.position="bottom", aspect.ratio=length(unique(Idents(sc_obj)))/length(features_extra_scn) + 0.1)


print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()