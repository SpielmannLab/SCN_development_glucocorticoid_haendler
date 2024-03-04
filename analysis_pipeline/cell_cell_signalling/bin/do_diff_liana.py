#!/usr/bin/env python
# Use Liana version 1.0.4. Installed in conda env liana_dev
# ----- Load libraries ------#
import argparse
import os

import decoupler as dc
import liana as li
import pandas as pd
import plotnine as p9
import scanpy as sc
# Set figure window
from matplotlib import rcParams
# Import DESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Define functions


def get_args():
    parser = argparse.ArgumentParser(description="A script to perform differential ")
    # Add arguments
    parser.add_argument(
        "--input_file",
        type=str,
        help="Path to the input file containing AnnData object with both the conditions to compare.",
    )
    parser.add_argument(
        "--celltype_column",
        type=str,
        help="The metadata column that has celltype information",
    )
    parser.add_argument(
        "--sample_column",
        type=str,
        help="The metadata column that has sample information. Usually orig.ident",
    )
    parser.add_argument(
        "--condition_column",
        type=str,
        help="The metadata column that has the information about the two conditions to be compared. Usually batch, treatment, replicate etc.",
    )
    parser.add_argument(
        "--condition_ref_level",
        type=str,
        help="Value of the reference condition. Is a value in the condition_column column",
    )
    parser.add_argument(
        "--condition_alt_level",
        type=str,
        help="Value of the alternate condition. Is a value in the condition_column column",
    )
    parser.add_argument(
        "--resource_name",
        default="mouseconsensus",
        type=str,
        help="LR resource to use for scoring. Run with --provide_more_help to figure out what you can use. Check logfile",
    )
    parser.add_argument(
        "--tileplot_width",
        default=12,
        type=int,
        help="Choose the width for tileplot. Default is 12",
    )
    parser.add_argument(
        "--tileplot_height",
        default=12,
        type=int,
        help="Choose the height for tileplot. Default is 12",
    )
    # Parse the command-line arguments and return
    args = parser.parse_args()
    return args


# Run DE analysis for each cluster at a time


def perform_de_on_pseudo_bulk(
    pdata, celltype_column, condition_column, condition_ref_level, condition_alt_level
):
    dea_results = {}
    quiet = True  # to supress verbosity

    for cell_group in pdata.obs[celltype_column].unique():
        # Select cell profiles
        ctdata = pdata[pdata.obs[celltype_column] == cell_group].copy()

        # Obtain genes that pass the edgeR-like thresholds
        # NOTE: QC thresholds might differ between cell types, consider applying them by cell type
        genes = dc.filter_by_expr(
            ctdata,
            group=condition_column,
            min_count=5,  # a minimum number of counts in a number of samples
            min_total_count=10,  # a minimum total number of reads across samples
        )

        # Filter by these genes
        ctdata = ctdata[:, genes].copy()

        # Build DESeq2 object
        # NOTE: this data is actually paired, so one could consider fitting the patient label as a confounder
        dds = DeseqDataSet(
            adata=ctdata,
            design_factors=condition_column,
            # set control as reference
            ref_level=[condition_column, condition_ref_level],
            refit_cooks=True,
            quiet=quiet,
        )

        # Compute LFCs
        dds.deseq2()
        # Contrast between condition_alt_level and condition_ref_level
        stat_res = DeseqStats(
            dds,
            contrast=[condition_column, condition_alt_level, condition_ref_level],
            quiet=quiet,
        )
        stat_res.quiet = quiet
        # Compute Wald test and display summary if needed
        stat_res.summary()
        # Shrink LFCs
        # {condition_column}_cond_vs_ref
        stat_res.lfc_shrink(
            coeff=f"{condition_column}_{condition_alt_level}_vs_{condition_ref_level}"
        )

        dea_results[cell_group] = stat_res.results_df

    return dea_results


def main():
    # ----- Get arguments from command line call ------#
    args = get_args()
    input_file = args.input_file
    celltype_column = args.celltype_column
    sample_column = args.sample_column
    condition_column = args.condition_column
    condition_ref_level = args.condition_ref_level
    condition_alt_level = args.condition_alt_level
    resource_name = args.resource_name
    tileplot_width = args.tileplot_width
    tileplot_height = args.tileplot_height

    # ----- Load data ------#
    # load a h5ad file
    adata = sc.read_h5ad(input_file)
    # first save the count matrix into a layer, because we will use raw counts to generate pseudobulk profiles
    # Note!!! adata.X should have normalized counts and adata.raw.X should have raw counts
    adata.layers["counts"] = adata.raw.X.copy()

    # Prepare folders for saving output
    os.makedirs(name="figures", exist_ok=True)
    os.makedirs(name="scores", exist_ok=True)

    # Generate UMAP for visualization
    rcParams["figure.figsize"] = (5, 5)
    sc.pl.umap(
        adata,
        color=[sample_column, celltype_column, condition_column],
        frameon=False,
        save="_sc_metadata",
        show=False,
        ncols=1,
    )

    # Generate pseudo-bulk data for DESeq2 analysis and plot how the pseudobulk data varies
    pdata = dc.get_pseudobulk(
        adata,
        sample_col=sample_column,
        groups_col=celltype_column,
        layer="counts",
        mode="sum",
        min_cells=10,
        min_counts=10000,
    )
    # Plot quality control to make sure that the pseudobulk profiles are not completely off between different metadata values
    dc.plot_psbulk_samples(
        pdata,
        groupby=[sample_column, celltype_column, condition_column],
        figsize=(20, 4),
        return_fig=False,
        save="figures/pseudo_bulk_qc.png",
    )

    # Run DE analysis for each cluster at a time
    dea_res = perform_de_on_pseudo_bulk(
        pdata,
        celltype_column,
        condition_column,
        condition_ref_level,
        condition_alt_level,
    )
    # Now concatenate the dataframe to contain DE info for all cell types together
    dea_df = pd.concat(dea_res)
    dea_df = (
        dea_df.reset_index()
        .rename(columns={"level_0": celltype_column})
        .rename(columns={"level_1": "index"})
        .set_index("index")
    )
    dea_df.head()

    # Now get the LR strenths for the Differentially Expressed ligand and receptors
    # They expression of ligands and receptors should change in the same direction.
    # That is both decrease or increase

    # Keep one condition - The alternate
    adata_alt = adata[adata.obs[condition_column] == condition_alt_level].copy()
    lr_res_alt = li.multi.df_to_lr(
        adata_alt,
        dea_df=dea_df,
        resource_name=resource_name,
        expr_prop=0.1,  # calculated for adata as passed - used to filter interactions
        groupby=celltype_column,
        stat_keys=["stat", "pvalue", "padj"],
        use_raw=False,
        complex_col="stat",  # NOTE: we use the Wald Stat to deal with complexes
        verbose=True,
        return_all_lrs=False,
    )
    # Save the unfiltered results to a file
    lr_res_alt.to_csv("scores/Dysregulated_LRs_unfiltered.tsv", sep="\t")
    # Now keep only the ones ones where the ligand, AND the receptor, AND the interaction change significantly
    lr_res_alt_significant = lr_res_alt[
        (lr_res_alt["ligand_padj"] < 0.05)
        & (lr_res_alt["receptor_padj"] < 0.05)
        & (lr_res_alt["interaction_padj"] < 0.05)
    ]
    # Save the filtered results to a file
    lr_res_alt_significant.to_csv(
        "scores/Dysregulated_LRs_only_significant_used_to_plot.tsv", sep="\t"
    )
    # Make a tileplot with top 40 LR pairs
    lr_tileplot = li.pl.tileplot(
        liana_res=lr_res_alt_significant,
        fill="stat",
        label="stat",
        label_fun=lambda x: round(x, ndigits=2),
        top_n=40,
        orderby="interaction_padj",
        orderby_ascending=True,
        orderby_absolute=False,
        source_title="Ligand",
        target_title="Receptor",
    )
    lr_tileplot = (
        lr_tileplot
        + p9.scale_fill_gradient2(
            low="#2166ac", high="#b2182b", mid="#f7f7f7", midpoint=0, limits=(-10, 10)
        )
        + p9.ggtitle(
            f"Wald statistic of change in expression in {condition_alt_level} compared to {condition_ref_level}"
        )
    )
    lr_tileplot.save(
        filename="figures/LR_TilePlot_significant_all_senders.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )
    # Now make a tile plot for every source cell type individually
    for source_cell in lr_res_alt_significant["source"].unique():
        lr_res_to_plot = lr_res_alt_significant[
            lr_res_alt_significant["source"] == source_cell
        ]
        lr_tileplot = li.pl.tileplot(
            liana_res=lr_res_to_plot,
            fill="stat",
            label="stat",
            label_fun=lambda x: round(x, ndigits=2),
            top_n=50,
            orderby="interaction_padj",
            orderby_ascending=True,
            orderby_absolute=False,
            source_title="Ligand",
            target_title="Receptor",
        )
        lr_tileplot = (
            lr_tileplot
            + p9.scale_fill_gradient2(
                low="#2166ac",
                high="#b2182b",
                mid="#f7f7f7",
                midpoint=0,
                limits=(-10, 10),
            )
            + p9.ggtitle(
                f"Wald statistic of change in expression in {condition_alt_level} compared to {condition_ref_level}"
            )
        )
        lr_tileplot.save(
            filename=f"figures/LR_TilePlot_significant_sender_{source_cell}.pdf",
            width=tileplot_width,
            height=tileplot_height,
        )


if __name__ == "__main__":
    main()
