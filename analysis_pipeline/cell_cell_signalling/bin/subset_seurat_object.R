#!/usr/bin/env Rscript
# Call from conda-environment MegaBundle
"
Usage: subset_seurat_object.R --file_sc_obj=<file> --subset_key=<value> --subset_values=<value>
Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		Seurat Object to be downsampled
  --subset_key=<value>    Metadata column name as key for subsetting
  --subset_values=<value> Values in the metadata column to subset (separate by <,> for multiple)
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))

# --- Parameters: Read
arguments <- docopt(doc, quoted_args = TRUE)
file_sc_obj <- arguments$file_sc_obj
subset_key <- arguments$subset_key
subset_values <- arguments$subset_values %>%
    strsplit(split=",") %>%
    unlist()

message("file_sc_obj: ", file_sc_obj)
message("subset_key: ", subset_key)
message("subset_values: ", paste(subset_values, collapse = ", "))

## Read the Seurat object
sc_obj <- readRDS(file_sc_obj)

print(sc_obj)

# Do subset only if subset_key found in the sc_obj metadata.
# Else, just output the distinct rows of the metadata
if ((subset_key %in% colnames(sc_obj@meta.data)) && !any(subset_values %in% c("null", "NULL")) && !any(is.null(subset_values))) {

    # Downsample using ths subset function in SeuratObject
    Idents(sc_obj) <- sc_obj[[subset_key]]
    sc_obj_subset <- subset(sc_obj, idents = subset_values)

    filename=gsub(file_sc_obj, pattern=".rds", replacement = paste0("_", subset_key, "_subsetted_.rds"))
    saveRDS(sc_obj_subset, file = filename)

} else if (!(subset_key %in% colnames(sc_obj@meta.data))) {

    table <- sc_obj@meta.data %>%
        head(n = 30)
    filename <- gsub(basename(file_sc_obj),
        pattern = ".rds",
        replacement = "_available_metadata.tsv")
    write.table(table, file = filename, sep = "\t", row.names = FALSE)
    warning("Check subset_key. A metadata table has been saved")

} else if ((subset_values %in% c("null", "NULL")) || is.null(subset_values)) {

    table <- sc_obj@meta.data %>%
        distinct(across(any_of(subset_key)), .keep_all = TRUE)
    filename <- gsub(basename(file_sc_obj),
        pattern = ".rds",
        replacement = "_available_metadata.tsv")
    write.table(table, file = filename, sep = "\t", row.names = FALSE)
    warning("Check subset_values. A metadata table has been saved")
}
