#!/usr/bin/env Rscript
# Call from conda-environment ggplot_essentials
" Do statistics on ssGSEA scores calculated by ESCAPE package - The output of this script includes include wilcoxon test and median values.
Usage: statistics_gsea.R --file_sc_gsea_metadata=<file> --metadata_to_analyse=<value> --groupby=<value> --groups_to_compare=<value> --splitby=<value> --multiple_comparison_adjustment_method=<value>

Options:
  -h --help			Show this screen.
  --file_sc_gsea_metadata=<file>		RDS file containing sc_obj metadata and ssGSEA scores
  --metadata_to_analyse=<value>  Use <ESCAPE_>. This has been appendend to all ssGSEA enrichment score colnames
  --groupby=<value>   Which metadata column of Seurat object to group for statistical testing
  --splitby=<value>     The statistical testing will be done on these separately
  --multiple_comparison_adjustment_method=<value>       Which algorimth to use for multiple testing adjustment
  --groups_to_compare=<value>   Which groups to compare pariwise: <WT_vs_Del,Dup_vs_Del,...>
" -> doc

suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(tidyr))
suppressPackageStartupMessages(library(rstatix))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

list2env(x = arguments, envir = environment())

# print out messages as well as deal with NULL values
message("file_sc_gsea_metadata: ", file_sc_gsea_metadata)
message("metadata_to_analyse: ", metadata_to_analyse)
message("groupby: ", groupby)
if (splitby %in% c("null", "NULL")) splitby <- NULL
message("splitby: ", splitby)
if (groups_to_compare %in% c("null", "NULL")) groups_to_compare <- NULL
if (!is.null(groups_to_compare)){
    groups_to_compare <- groups_to_compare %>% 
        strsplit(split=",") %>% 
        unlist() %>% 
        lapply(., FUN = function(x){strsplit(x, split = "_vs_")}) %>% 
        unlist(recursive = FALSE)
}
message("groups_to_compare: ", paste(groups_to_compare, collapse = ", "))
message("multiple_comparison_adjustment_method: ", multiple_comparison_adjustment_method)

# --- Read files
message("Loading Data")
sc_gsea_metadata <- readRDS(file_sc_gsea_metadata)

# groupby <- "condition"
# splitby <- "integrated_snn_res.0.05"
# metadata_to_analyse <- "ESCAPE_"

enrichments <- cbind(select(sc_gsea_metadata, all_of(c(groupby, splitby))),
    select(sc_gsea_metadata, starts_with(metadata_to_analyse)))
enrichments <- pivot_longer(enrichments, cols = matches(metadata_to_analyse), names_to = "geneset", values_to = "score")

# Change the name of the groupby column so it is easy to manipulate.
colnames(enrichments)[colnames(enrichments) == groupby] <- "group"

message("Running statistics")
# Do the stats and median calculations based on whether splitby is provided or not.
if(is.null(splitby)) {
    # Do statistics
    stat.test <- enrichments %>%
        group_by(geneset) %>%
        wilcox_test(score ~ group, comparisons = groups_to_compare, p.adjust.method = multiple_comparison_adjustment_method) %>%
        adjust_pvalue(method = multiple_comparison_adjustment_method)  %>%
        add_significance()
    # Calculate Medians
    medians <-  enrichments %>%
        group_by(geneset, group) %>%
        summarize(median_ssGSEA_score = median(score, na.rm = TRUE))
    # Combine statistics with median ssGSEA score
    stat.test_with_median <- left_join(x = stat.test, y = medians, by = join_by(geneset == geneset, group1 == group)) %>%
        left_join(y = medians, by = join_by(geneset == geneset, group2 == group), suffix=c("1","2"))
} else {
    colnames(enrichments)[colnames(enrichments) == splitby] <- "split"
    # Do statistics
    stat.test <- enrichments %>%
        group_by(split, geneset) %>%
        wilcox_test(score ~ group, comparisons = groups_to_compare, p.adjust.method = multiple_comparison_adjustment_method) %>%
        adjust_pvalue(method = multiple_comparison_adjustment_method)  %>%
        add_significance()
    # Calculate Medians
    medians <-  enrichments %>%
        group_by(split, geneset, group) %>%
        summarize(median_ssGSEA_score = median(score, na.rm = TRUE))
    # Combine statistics with median ssGSEA score
    stat.test_with_median <- left_join(x = stat.test, y = medians, by = join_by(split == split, geneset == geneset, group1 == group)) %>%
        left_join(y = medians, by = join_by(split == split, geneset == geneset, group2 == group), suffix=c("1","2"))
}
# Save files
filename <- paste0("wilcox_test_statistics_with_median_",
    gsub(basename(file_sc_gsea_metadata),
        pattern = ".rds",
        replacement = paste0("_", groupby, "_", splitby, "_.tsv")))

message("Saving stat results")
write.table(stat.test_with_median, file = filename, col.names = TRUE, sep = "\t", row.names = FALSE)