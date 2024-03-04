#!/usr/bin/env Rscript
# Call from conda-environment 
" Make Violin plots based on the GSEA analysis results

Usage: vlnplot_gsea_per_pathway.R --file_sc_gsea_metadata=<file> --metadata_to_plot=<value> --vln_groupby=<value> --vln_splitby=<value> --pt_size=<value> --vln_width=<value> --vln_height=<value> --vln_colors=<value>

Options:
    -h --help			Show this screen.
    --file_sc_gsea_metadata=<file>		RDS file containing sc_obj metadata and ssGSEA scores
    --metadata_to_plot=<value>  Use <ESCAPE_>. This has been appendend to all ssGSEA enrichment score colnames
    --vln_groupby=<value>   Which metadata column to use as the x-axis in vlnplot
    --vln_splitby=<value>   Which metadata column to use for coloring the vlnplot and for stat testing.
    --pt_size=<value>       Size for the points on the violin plot
    --vln_width=<value>     Width of the plot output
    --vln_height=<value>    Height of the plot output
    --vln_colors=<value>    Colors for the violin plot
" -> doc

# --- load libraries
suppressMessages(library(docopt))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrastr))


# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

list2env(x = arguments, envir = environment())

pt_size <- as.numeric(pt_size)
vln_width <- as.numeric(vln_width)
vln_height <- as.numeric(vln_height)
vln_colors <- vln_colors %>% 
    strsplit(split = ",") %>%
    unlist()

message("file_sc_gsea_metadata: ", file_sc_gsea_metadata)
message("metadata_to_plot: ", metadata_to_plot)
message("vln_groupby: ", vln_groupby)
message("vln_splitby: ", vln_splitby)
message("pt_size: ", pt_size)
message("vln_width: ", vln_width)
message("vln_height: ", vln_height)
message("vln_colors: ", vln_colors)


# --- Read files
message("Loading Data")
sc_gsea_metadata <- readRDS(file_sc_gsea_metadata)

pathways <- grep(colnames(sc_gsea_metadata), pattern = metadata_to_plot, value = TRUE)
for (pathway in pathways) {
    sc_gsea_metadata_subset <- sc_gsea_metadata %>%
        select(all_of(c(vln_groupby, vln_splitby, pathway)))
    names(sc_gsea_metadata_subset) <- c("groupby", "splitby", "ssGSEA_score")

    ## Remove groups and splits which have fewer than 2 datapoints
    sc_gsea_metadata_subset <- sc_gsea_metadata_subset %>%
        group_by(groupby, splitby) %>%
        filter(n() >= 2)

    ## Remove all groups in which there are more than one unique value in splitby (to do pairwise statistical testing)
    groups_without_multiple_splitbys <- sc_gsea_metadata_subset %>% 
        group_by(groupby) %>%
        summarize(n_unique = length(unique(splitby))) %>%
        filter(n_unique ==1) %>% pull(groupby)
    sc_gsea_metadata_subset <- sc_gsea_metadata_subset %>% 
        filter(!(groupby) %in% groups_without_multiple_splitbys)


    # Do the Wilcox test
    stat_test_results <- sc_gsea_metadata_subset %>% 
        group_by(groupby) %>%
        wilcox_test(ssGSEA_score ~ splitby, p.adjust.method = "bonferroni")
    filename = paste0("wilcox_test_statistics_",
        gsub(basename(file_sc_gsea_metadata),
            pattern=".rds",
            replacement=paste0("_", vln_groupby, "_", vln_splitby, "_", pathway, ".tsv")))
    write.table(stat_test_results, file = filename, col.names = TRUE, sep = "\t", row.names = FALSE)


    vln_plot <- ggplot(sc_gsea_metadata_subset, aes(y=ssGSEA_score, x=groupby, color = splitby)) +
        geom_violin(position = position_dodge(width = 0.9), draw_quantiles = 0.5) +
        geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), size = pt_size, shape = 16) +
        scale_fill_manual(values = vln_colors) +
        scale_color_manual(values = vln_colors) +
            theme_classic() +
            theme(strip.background=element_blank(),
                text=element_text(size=10, family="sans"),
                legend.position="bottom",
                axis.text=element_text(size=10),
                axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
                strip.text=element_text(size = 8)) +
            xlab("") +
            ggtitle(pathway)
    
    # Whether to work with pvalue or pvalue-adjusted (which is only present in case of multiple comparisons)?
    if ("p.adj" %in% colnames(stat_test_results)) {
        label_on_plot <- "p.adj"
    } else {label_on_plot <- "p"}

    # Add the caluclated p-value to the vln_plot
    stat_test_results <- stat_test_results %>%
        filter(.data[[label_on_plot]] < 0.05) %>%
        add_xy_position(x = "groupby", dodge = 0.9)
    
    vln_plot <- vln_plot + stat_pvalue_manual(stat_test_results, label = label_on_plot)

    filename = paste0("Violinplot_",
        gsub(basename(file_sc_gsea_metadata),
            pattern=".rds",
            replacement=paste0("_", vln_groupby, "_", vln_splitby, "_", pathway, ".pdf")))
    ggsave(vln_plot, filename=filename, width = vln_width, height = vln_height)
}