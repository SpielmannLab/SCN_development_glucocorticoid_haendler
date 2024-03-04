#!/usr/bin/env python
# Use Liana version 1.0.4. Installed in conda env liana_dev
# ----- Load libraries ------#
import argparse
import os

import liana as li
import pandas as pd
import plotnine as p9
import numpy as np

# Define parameters
file_diff_liana_scores_filtered = "/Users/sreenivasan/Documents/Works/AstizSCN/cell_cell_signalling_analysis/cell_signalling_output/liana_custom_lr/sig_lr_in_collecTRI_Nr3c1_regulon.tsv"
lr_res_alt_filtered = pd.read_csv(file_diff_liana_scores_filtered, sep="\t")
unique_id = "collectTRI_Nr3c1_regulons_super_filtered"

os.makedirs(name="figures", exist_ok=True)
# Create a column with ordering based on the biggest change in either receptor
# or ligand expression
lr_res_alt_filtered["ordering"] = lr_res_alt_filtered.ligand_stat.combine(
    lr_res_alt_filtered.receptor_stat, lambda x, y: max(abs(x), abs(y))
)

# Make a tileplot with top 40 LR pairs
lr_tileplot = li.pl.tileplot(
    liana_res=lr_res_alt_filtered,
    fill="stat",
    label="padj",
    label_fun = lambda x: '*' if x < 0.05 else np.nan,
    top_n=40,
    orderby="ordering",
    orderby_ascending=False,
    orderby_absolute=False,
    source_title="Ligand",
    target_title="Receptor",
)
lr_tileplot = (
    lr_tileplot
    + p9.scale_fill_gradient2(
        low="#2166ac", high="#b2182b", mid="#f7f7f7", midpoint=0, limits=(-15, 15)
    )
    + p9.ggtitle("Wald statistic of change in expression")
    + p9.theme_classic()
    + p9.theme(axis_text_x = p9.element_text(rotation = 90))
)

tileplot_width = 5
tileplot_height = 5
lr_tileplot.save(
    filename=f"figures/LR_TilePlot_filtered_{unique_id}.pdf",
    width=tileplot_width,
    height=tileplot_height,
)
