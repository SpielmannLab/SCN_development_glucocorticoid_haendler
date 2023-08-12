#!/usr/bin/env Rscript
" Transfer the pseudotime and UMAP co-ordinates into the original Seurat object and plot here.
Usage: pseudotime_in_seurat.R --file_sc_obj=<file> --filename_prefix=<value> --assay=<value> --file_cds_obj=<file> --genes=<value> --gex_pt_size=<value> --gex_width=<value> --gex_height=<value>

Options:
    -h --help			Show this screen.
    --file_sc_obj=<file>		The rds file of a seurat object to convert to monocle object
    --filename_prefix=<value>		Optional. Provide a string to prefix all filenames. If <NA>, uses unique value from annotations
    --assay=<value>                 Which assay to use 
    --file_cds_obj=<file>		    The rds file of a monocle object (CDS) to perform trajectory analysis
    --genes=<value>                 Genes to plot the expression values as a function of pseudotime
    --gex_pt_size=<value>           Size of dots in the gex vs pseudotime plot
    --gex_width=<value>             Size of the gex vs pseudotime plot
    --gex_height=<value>            Size of the gex vs pseudotime plot
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

list2env(x = arguments, envir = environment())
genes <- genes %>%
    strsplit(split = ",") %>% 
    unlist()
gex_pt_size <- as.numeric(gex_pt_size)
gex_width <- as.numeric(gex_width)
gex_height <- as.numeric(gex_height)

if ((filename_prefix %in% c("null", "NULL") || is.null(filename_prefix))) {
    filename_prefix <- gsub(file_sc_obj, pattern = ".rds", replacement = "")
}

message("file_sc_obj: ", file_sc_obj)
message("assay: ", assay)
message("file_cds_obj: ", file_cds_obj)
message("genes: ", paste(genes, collapse = ", "))
message("gex_pt_size: ", gex_pt_size)
message("gex_width: ", gex_width)
message("gex_height: ", gex_height)

# read file and apply parameters
sc_obj <- readRDS(file_sc_obj)
DefaultAssay(sc_obj) <- assay

cds <- readRDS(file_cds_obj)

# Extract the pseudotime information from cds and add it to Seurat object as new Metadata
sc_obj_subset <- subset(sc_obj, cells = Cells(cds))
if(!identical(Cells(sc_obj_subset), Cells(cds))) stop("The cells in cds obj does not match with that in the seurat object for pseudotime transfer")

sc_obj_subset$pseudotime <- pseudotime(cds)

filename <- paste0(filename_prefix, "_traj_gex_in_pseudotime.rds")
saveRDS(sc_obj_subset, file = filename)


# Make plot
group_by <- "developmental_timepoint"
colors <- c("#a1dab4","#41b6c4","#2c7fb8","#253494")


for (gene in genes) {
    plot <- FeatureScatter(sc_obj_subset,
            feature1 = "pseudotime",
            feature2 = gene,
            slot = "data",
            pt.size = gex_pt_size,
            jitter = FALSE,
            group.by = group_by,
            cols = colors,
            plot.cor = FALSE) +
        geom_smooth(formula = y ~ splines::ns(x, df=5),
            se = FALSE,
            color = "black") +
        ylim(0.5, NA) +
        theme_classic() +
        theme(aspect.ratio = 0.5)

    filename <- paste0(filename_prefix, "_traj_gex_in_pseudotime_", gene, ".pdf")
    ggsave(plot = plot, filename = filename, width = gex_width, height = gex_height)
}