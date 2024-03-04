#!/usr/bin/env Rscript
# Call ggplot_essentials conda environment
" Make volcano plot from statistics ouput of ssGSEA scores.

Usage: volcano_gsea.R --file_gsea_statistics_with_median=<file> --width=<value> --height=<value> --pt_size=<value>

Options:
    -h --help			Show this screen.
    --file_gsea_statistics_with_median=<file>		A tsv file containing the GSEA stats
    --width=<value>                                 Width of the plot
    --height=<value>                                Height of the plot 
    --pt_size=<value>                               Size of points
" -> doc

# --- load libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

list2env(x = arguments, envir = environment())

width <- as.numeric(width)
height <- as.numeric(height)
pt_size <- as.numeric(pt_size)

message("file_gsea_statistics_with_median: ", file_gsea_statistics_with_median)
message("width: ", width)
message("height: ", height)
message("pt_size: ", pt_size)

# --- Read files
message("Loading Data")
stat.test_with_median <- read.table(file_gsea_statistics_with_median, sep="\t", header = TRUE)
# genesetColorsLabels <- read.table(file_genesetColorsLabels, sep = "\t", header = TRUE)

stat.test_with_median <- stat.test_with_median %>% 
    mutate(diff_medians = (median_ssGSEA_score1 - median_ssGSEA_score2)) %>%
    mutate(comparison = paste0(group1, "_vs_", group2))

filename=paste0("Volcano_",gsub(basename(file_gsea_statistics_with_median), pattern=".tsv", replacement=".pdf"))
    plot <- ggplot(data = stat.test_with_median,
        aes(x = diff_medians, y = -log10(p.adj))) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_point(size = pt_size) +
    theme_classic() +
    xlab("Change in median ssGSEA score group1-group2") +
    theme(aspect.ratio=1,
        axis.title=element_text(size=10),
        axis.text=element_text(size=10))

if ("split" %in% colnames(stat.test_with_median)) {
    plot <- plot + facet_grid(comparison ~ split, scales = "free")
} else {
    plot <- plot + facet_grid(comparison ~ ., scales = "free")
}

ggsave(plot = plot, filename = filename, width = width, height = height)
