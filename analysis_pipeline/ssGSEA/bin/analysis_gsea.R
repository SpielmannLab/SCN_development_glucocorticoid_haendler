#!/usr/bin/env Rscript
# Call from conda-environment ssGSEA
" Perform a GSEA analysis on the submitted Seurat object, and output the Seurat Object with enrichment scores.
It will not run again, if the file already exists in the outfolder

Usage: runGSEA.R --file_sc_obj=<file> --geneSet=<value> --assay=<value> --species=<value> --outfolder=<path>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		Desired Seurat Object
  --geneSet=<value>  Library,NamePattern Check www.gsea-msigdb.org/gsea/msigdb/collections.jsp
  --assay=<value>   Assay in Seurat object to use <RNA>, <spliced>, or <unspliced>
  --species=<value>   Homo sapiens or Mus musculus
  --outfolder=<value> To check if the rds file already exists in the output folder
" -> doc

# --- load libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(parallel))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args = TRUE)
print(arguments)

list2env(x = arguments, envir = environment())

geneSet <- strsplit(geneSet, split = ",") %>%
  unlist()
geneSetLibrary <- geneSet[1]
geneSetNamePattern <- geneSet[2]

message("file_sc_obj: ", file_sc_obj)
message("geneSetLibrary: ", geneSetLibrary)
message("geneSetNamePattern: ", geneSetNamePattern)
message("assay: ", assay)
message("species: ", species)
message("outfolder: ", outfolder)

# parameter set
ncores <- parallel::detectCores()
message("ncores: ", ncores)

# Define functions
getGeneSets_from_msigdb <- function(geneSetLibrary, geneSetNamePattern, species) {
  if (is.na(geneSetNamePattern)) {
    selected_gene_sets <- NULL
  } else {
    # -- Pull gene sets from msigdbr
    all_gene_sets <- msigdbr(species = species) # Get all gene sets in the MSIG database
    selected_gene_sets <- all_gene_sets %>%
      filter(gs_cat == geneSetLibrary) %>%
      filter(grepl(gs_name, pattern = geneSetNamePattern, ignore.case = TRUE)) %>%
      distinct(gs_name) %>%
      .$gs_name
  }
  gene.sets <- getGeneSets(library = geneSetLibrary, gene.sets = selected_gene_sets, species = species) # Get the user defined gene sets
}

getGeneSets_from_file <- function(geneSetFile) {
  lines <- readLines(geneSetFile) %>%
    strsplit(split = ":")
  gene.sets <- lapply(X = lines, FUN = function(line) {
    set <- strsplit(line[2], split = ",") %>%
      unlist() %>%
      list()
    names(set) <- line[1]
    return(set)
  })
  gene.sets <- unlist(gene.sets, recursive = FALSE)
  return(gene.sets)
}

# If the provide geneSet is a file, load it, else, pull from msigdb
if (geneSetLibrary == "custom") {
  # If a custom file is provided, 
  # Set the namepattern to something meaningful
  # And set the filename to something meaningful
  geneSetFile <- geneSetNamePattern
  geneSetNamePattern <- "file"
  gene.sets <- getGeneSets_from_file(geneSetFile = geneSetFile)
} else {
  gene.sets <- getGeneSets_from_msigdb(
    geneSetLibrary = geneSetLibrary,
    geneSetNamePattern = geneSetNamePattern,
    species = species)
}

# Check if the the gsea for this geneset has already been done. If so, do not repeat
filename <- gsub(file_sc_obj, pattern = ".rds", replacement = paste0("_", assay, "_", geneSetLibrary, "_", geneSetNamePattern, ".rds"))
file_exists <- list.files(path = outfolder, pattern = filename) %>%
  length() %>%
  as.logical()

if (file_exists) {
  warning("This GSEA analysis has already been done. So, aborting the call")
  file.copy(from = paste0(outfolder, "/", filename), to = filename)
  quit("no")
}

# --- Read files
sc_obj <- readRDS(file_sc_obj)
DefaultAssay(sc_obj) <- assay

sc_obj <- sc_obj %>%
  NormalizeData(normalization.method = "LogNormalize")



# calculate enrichment scores for the cells using ESCAPE
start_time <- Sys.time()
ES <- enrichIt(
  obj = sc_obj,
  gene.sets = gene.sets,
  groups = 10000, cores = ncores
)
message("Time for ssGSEA calculations: ", Sys.time() - start_time)
colnames(ES) <- paste0("ESCAPE_", colnames(ES)) # This is done so that these columns can be selected in the plotting R-script

# Add the enrichment scores as metadata to the original sc_obj
sc_obj <- AddMetaData(sc_obj, ES)
sc_gsea_metadata <- sc_obj@meta.data

# save only the metadata
filename <- gsub(file_sc_obj, pattern = ".rds", replacement = paste0("_", assay, "_", geneSetLibrary, "_", geneSetNamePattern, ".rds"))
saveRDS(sc_gsea_metadata, file = filename)
