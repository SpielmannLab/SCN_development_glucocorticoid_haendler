import scvelo as scv
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import sys

# Get arguments
# path="/data/humangen_mouse/mmca/sox9_reg_ko"
#file_h5ad="/data/humangen_mouse/mmca/sox9_reg_ko/obj_sox9_wt_WT.h5ad"
#plot_filetype=".pdf"
#proportions_groupby="sub_trajectory"

file_h5ad=sys.argv[1]
plot_filetype=sys.argv[2]
proportions_groupby=sys.argv[3]

print("file_h5ad="+file_h5ad)
print("plot_filetype="+plot_filetype)
print("proportions_groupby="+proportions_groupby)

adata = scv.read(file_h5ad, cache=True) # This is the file created by seurat_to_anndata.R code
# Make sure that the counts_matrix going to the adata$X is the "Data" slot from the RNA assay
resolution =[meta for meta in adata.obs_keys() if proportions_groupby in meta]
# QC-ensure the importw worked well
print(type(resolution[0]))

plot_filename=os.path.basename(file_h5ad)
print(plot_filename)
plot_filename="_" + os.path.splitext(plot_filename)[0] + plot_filetype
print(plot_filename)
sc.pl.umap(adata, color=resolution[0], frameon=False, save=plot_filename)


#
# Normalize and PCA
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=50, n_neighbors=30)

# Calculate velocity and save plot
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

plot_filename=os.path.basename(file_h5ad)
plot_filename=os.path.splitext(plot_filename)[0] + plot_filetype
scv.pl.velocity_embedding_stream(adata, basis='umap', color=resolution[0], legend_loc="right margin", save=plot_filename)

# QC - Plot proportions of spliced and unspliced transcripts
plot_filename=os.path.basename(file_h5ad)
plot_filename=os.path.splitext(plot_filename)[0] + plot_filetype
scv.pl.proportions(adata, groupby=str(resolution[0]), save=plot_filename)
