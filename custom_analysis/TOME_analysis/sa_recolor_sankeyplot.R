# conda activate AstizSCN

library(Seurat)
gene <- "Nr3c1"
assay <- "RNA"
files <- list.files(pattern = ".rds", path = "/data/humangen_mouse/2023_glucocorticoid/output/TOMEbycell_state/processed")
annotations_column <- "Anno" 

# Get the expression of the genes for all the Seurat objects.
expression_list <- lapply(files, FUN = function(file_sc_obj){
    sc_obj <- readRDS(file_sc_obj)
    expression <- AverageExpression(sc_obj, 
        features = gene,
        group.by = annotations_column,
        assays = assay)
    expression <- expression[[assay]] %>%
        t() %>%
        data.frame() %>%
        tibble::rownames_to_column(var = "annotations_column")
})
expression <- do.call(expression_list, what = rbind)
names(expression) <- c("annotations_column", "lognorm")
expression$zero_to_one <- expression$lognorm/max(expression$lognorm)

# Get color codes for these epxressions for a given color gradient
colfunc <- colorRamp(c("#fff7bc", "#fec44f", "#d95f0e"))
expression$cols <- colfunc(expression$zero_to_one)
colorcodes <- rgb(expression$cols[,1],expression$cols[,2],expression$cols[,3], alpha = 255, maxColorValue = 255) %>% paste(collapse = '", "') 
colorcodes <- paste0('"',colorcodes,'"')

annotations <- paste(expression$annotations_column, collapse = '", "')
annotations <- paste0('"', annotations, '"')

# Recoloring Sankey plot with the gradient obtained with expression
# Sankeyplot color examples at:
# https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html
# Manually assign colors obtained above by doing <cat colorcodes>
my_color <- 'd3.scaleOrdinal() .domain(["GD17_?", "GD17_AVP GABAergic neuron immature", "GD17_Enk GABAergic/dopaminergic immature neuron", "GD17_GABAergic immature neuron", "GD17_GABAergic progenitor", "GD17_Glutamatergic neuron immature", "GD17_Neuroblast (astrocyte fate)", "GD17_Neuroblast (OPCs fate)", "GD17_Prlr GABAergic immature neurons", "GD17_Progenitor Enk GABAergic neurons", "GD17_Progenitor GABAergic neurons", "GD17_Progenitor Prlr GABAergic neurons", "GD17_Retinal ganglia cells", "GD17_Retinal ganglion cells", "PND10_Astrocytes", "PND10_AVP GABAergic neurons", "PND10_AVP/Prok2/Nms GABAergic neurons", "PND10_Enk/Prlr GABAergic neuron", "PND10_Ependymal cells", "PND10_GABAergic neurons", "PND10_Mature Oligodendrocytes", "PND10_Microglia", "PND10_Neuron (not clear which one)", "PND10_Oligodendrocytes precursors/NG2 cells", "PND10_Prlr GABAergic immature neurons", "PND10_Prlr GABAergic neuron", "PND10_VIP/Grp neurons (GABAergic)", "PND2_?", "PND2_AVP GABAergic neurons", "PND2_Bergmann glia", "PND2_Enk immature GABAergic neurons", "PND2_Ependymal cells", "PND2_GABAergic immature neurons", "PND2_GABAergic neurons", "PND2_Immature astrocytes", "PND2_Immature astrocytes/Bregmann glia", "PND2_Immature GABAergic/dopaminergic neurons", "PND2_Immature interneurons", "PND2_Intermediate progenitor cells/Prlr GABAergic neuron", "PND2_Intermediate progenitors/GABAergic", "PND2_Microglia", "PND2_Neuroendocrine cells", "PND2_Neuroendocrine cells (immature GABAergic neuron)", "PND2_Oligodendrocytes/NG2 cells", "PND2_Prlr GABAergic immature neurons", "PND2_Prlr Immature GABAergic/dopaminergic neurons", "PND2_VIP GABAergic neurons", "PND30_Astrocytes", "PND30_AVP GABAergic neurons", "PND30_AVP neurons (GABAergic neuron)", "PND30_CCK neurons (GABAergic)", "PND30_CCK/dopaminergic  (GABAergic)", "PND30_Enk GABAergic neurons", "PND30_Erk GABAergic neurons", "PND30_Mature Oligodendrocytes", "PND30_Microglia", "PND30_Oligodendrocytes progenitor cells or NG2", "PND30_Prlr GABAergic neurons", "PND30_VIP/GRP neurons (GABAergic)", "PND30_VIP/Prok2/Nms GABAergic neurons"]) .range(["#FEF3B5FF", "#FEF3B4FF", "#F2A53BFF", "#FEDF88FF", "#FEEFACFF", "#FEE89DFF", "#FEE89CFF", "#FEE28FFF", "#FEE89CFF", "#FEF5B8FF", "#FEEEA9FF", "#FEE79AFF", "#FEF3B4FF", "#FEF4B6FF", "#FBBC4AFF", "#FEECA6FF", "#FEE99FFF", "#FEE89DFF", "#F09F37FF", "#D95F0EFF", "#FEDC82FF", "#E47D21FF", "#FEEDA8FF", "#FECD63FF", "#FED36FFF", "#FEE18DFF", "#FEF3B4FF", "#FEDB81FF", "#FEF0AEFF", "#FECB5FFF", "#FEE290FF", "#FED069FF", "#FEE697FF", "#FCBF4CFF", "#FEE18EFF", "#FED87BFF", "#FEEAA1FF", "#FED574FF", "#FEE392FF", "#FEDC82FF", "#F7B142FF", "#F6AF41FF", "#F09D36FF", "#FEDB82FF", "#FEE698FF", "#FEDC82FF", "#FEF2B3FF", "#FCBF4CFF", "#FEEDA8FF", "#FEE697FF", "#FEF0AEFF", "#FEF2B2FF", "#FEF2B1FF", "#FEDD86FF", "#FEC552FF", "#FDC24DFF", "#FED26DFF", "#FED97CFF", "#FEEFACFF", "#FEF2B2FF"])'

p <- sankeyNetwork(Links = combined_sankey_links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name",
                   colourScale = my_color,
                   sinksRight=FALSE, fontSize=20)

filename <- "test.html"
saveWidget(p, file=filename)