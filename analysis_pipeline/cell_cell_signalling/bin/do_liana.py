#!/usr/bin/env python
# Use Liana version 1.0.4. Installed in conda env liana_dev
# ----- Load libraries ------#
# for error logging and informational messages
# to parse arguments
import argparse
import logging
import os

# for liana analysis and individual methods
import liana as li
# for visualization and toy data
import scanpy as sc
from liana.method import (cellchat, cellphonedb, connectome, geometric_mean,
                          logfc, natmi, singlecellsignalr)

logging.basicConfig(
    level=logging.INFO, filename="liana.log", format="%(levelname)s:%(message)s"
)


def get_args():
    parser = argparse.ArgumentParser(
        description="A comprehensive script to run multiple Liana methods. Check the liana.log file after run completion. Requires liana version 1.0.4 or higher"
    )
    # Add arguments
    parser.add_argument(
        "--input_file",
        type=str,
        default="test_data",
        help="Path to the input file containing AnnData object. Leave empty to get trial results for pbmc68k_reduced dataset.",
    )
    parser.add_argument(
        "--provide_more_help",
        action="count",
        help="Provide me more Liana specific help",
    )
    parser.add_argument(
        "--celltype_column",
        type=str,
        default="",
        help="The metadata column that has celltype information",
    )
    parser.add_argument(
        "--source_celltypes",
        type=str,
        default=None,
        help="Optional. The specific cell types to be plotted as signalling source. Comma separated",
    )
    parser.add_argument(
        "--target_celltypes",
        type=str,
        default=None,
        help="Optional. The specific cell types to be plotted as signalling target. Comma separated",
    )
    parser.add_argument(
        "--methodsToRun",
        default="CellPhoneDB",
        type=str,
        help="Methods to run using Liana. Run with --provide_more_help to figure out what you can use. Check logfile",
    )
    parser.add_argument(
        "--resource_name",
        default="consensus",
        type=str,
        help="LR resource to use for scoring. Run with --provide_more_help to figure out what you can use. Check logfile",
    )
    parser.add_argument(
        "--de_method",
        default="t-test",
        type=str,
        help="Choose from t-test_overestim_var, wilcoxon, logreg, t-test (Default)",
    )
    parser.add_argument(
        "--dotplot_width",
        default=12,
        type=int,
        help="Choose the width for dotplot. Default is 12",
    )
    parser.add_argument(
        "--dotplot_height",
        default=18,
        type=int,
        help="Choose the height for dotplot. Default is 12",
    )
    parser.add_argument(
        "--tileplot_width",
        default=12,
        type=int,
        help="Choose the width for tileplot. Default is 12",
    )
    parser.add_argument(
        "--tileplot_height",
        default=18,
        type=int,
        help="Choose the height for tileplot. Default is 12",
    )
    # Parse the command-line arguments and return
    args = parser.parse_args()
    return args


def run_geometric_mean(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using Geometric Means method. Check liana.log for more information."
    )
    logging.info(geometric_mean.describe())
    # Define to use global adata
    geometric_mean(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        verbose=True,
        key_added="geom_res",
        de_method=de_method,
    )
    # write the GeometricMean output to file
    geom_res = adata.uns["geom_res"]
    geom_res.to_csv("scores/GeometricMean_LR_scores.tsv", sep="\t")
    # Make dotplot
    geom_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="lr_gmeans",
        size="gmean_pvals",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=True,
        top_n=40,
        orderby="gmean_pvals",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(12, 12),
        # we filter the pvals column to <= 0.05
        filter_fun=lambda x: x["gmean_pvals"] <= 0.05,
        # uns_key to use, default is "liana_res"
        uns_key="geom_res",
    )
    geom_lr_dotplot.save(
        filename="figures/GeometricMean_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    # Make tileplot
    geom_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="means",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        filter_fun=lambda x: x["gmean_pvals"] <= 0.05,
        top_n=40,
        orderby="gmean_pvals",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="geom_res",
        figure_size=(12, 12),
    )
    geom_lr_tileplot.save(
        filename="figures/GeometricMean_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_cellphonedb(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using CellPhoneDB. Check liana.log for more information."
    )
    logging.info(cellphonedb.describe())
    cellphonedb(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        verbose=True,
        key_added="cpdb_res",
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    cpdb_res = adata.uns["cpdb_res"]
    cpdb_res.to_csv("scores/CellPhoneDB_LR_scores.tsv", sep="\t")
    # Make dotplot
    cpdb_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="lr_means",
        size="cellphone_pvals",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=True,
        top_n=40,
        orderby="cellphone_pvals",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(12, 12),
        # finally, since cpdbv2 suggests using a filter to FPs
        # we filter the pvals column to <= 0.05
        filter_fun=lambda x: x["cellphone_pvals"] <= 0.05,
        # uns_key to use, default is "liana_res"
        uns_key="cpdb_res",
    )
    cpdb_lr_dotplot.save(
        filename="figures/CellPhoneDB_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    # Make tileplot
    cpdb_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="means",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        filter_fun=lambda x: x["cellphone_pvals"] <= 0.05,
        top_n=40,
        orderby="cellphone_pvals",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="cpdb_res",  # NOTE: default is 'liana_res'
        figure_size=(12, 12),
    )
    cpdb_lr_tileplot.save(
        filename="figures/CellPhoneDB_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_connectome(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using Connectome method. Check liana.log for more information."
    )
    logging.info(connectome.describe())
    connectome(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        verbose=True,
        key_added="cnctm_res",
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    cnctm_res = adata.uns["cnctm_res"]
    cnctm_res.to_csv("scores/Connectome_LR_scores.tsv", sep="\t")
    # Make dotplot
    cnctm_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="expr_prod",
        size="scaled_weight",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=False,
        top_n=40,
        orderby="scaled_weight",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(12, 12),
        # uns_key to use, default is "liana_res"
        uns_key="cnctm_res",
    )
    cnctm_lr_dotplot.save(
        filename="figures/Connectome_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    # Make tileplot
    cnctm_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="means",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        top_n=40,
        orderby="scaled_weight",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="cnctm_res",  # NOTE: default is 'liana_res'
        figure_size=(12, 12),
    )
    cnctm_lr_tileplot.save(
        filename="figures/Connectome_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_natmi(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using NATMI method. Check liana.log for more information."
    )
    logging.info(natmi.describe())
    natmi(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        base=2,
        verbose=True,
        key_added="natmi_res",
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    natmi_res = adata.uns["natmi_res"]
    natmi_res.to_csv("scores/NATMI_LR_scores.tsv", sep="\t")
    # Make dotplot
    natmi_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="expr_prod",
        size="spec_weight",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=False,
        top_n=40,
        orderby="spec_weight",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(12, 12),
        # uns_key to use, default is "liana_res"
        uns_key="natmi_res",
    )
    natmi_lr_dotplot.save(
        filename="figures/NATMI_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    natmi_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="means",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        top_n=40,
        orderby="spec_weight",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="natmi_res",
        figure_size=(12, 12),
    )
    natmi_lr_tileplot.save(
        filename="figures/NATMI_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_singlecellsignalr(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using SingleCellSignalR method. Check liana.log for more information."
    )
    logging.info(singlecellsignalr.describe())
    singlecellsignalr(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        base=2,
        verbose=True,
        key_added="sscr_res",
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    sscr_res = adata.uns["sscr_res"]
    sscr_res.to_csv("scores/SingleCellSignalR_LR_scores.tsv", sep="\t")
    # Make dotplot
    sscr_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="lrscore",
        size="lrscore",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=False,
        top_n=40,
        orderby="lrscore",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(8, 7),
        # uns_key to use, default is "liana_res"
        uns_key="sscr_res",
    )
    sscr_lr_dotplot.save(
        filename="figures/SingleCellSignalR_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    sscr_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="means",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        top_n=40,
        orderby="lrscore",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="sscr_res",
        figure_size=(12, 12),
    )
    sscr_lr_tileplot.save(
        filename="figures/SingleCellSignalR_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_logfc(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using log2FC method. Check liana.log for more information."
    )
    logging.info(logfc.describe())
    logfc(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        base=2,
        verbose=True,
        key_added="logfc_res",
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    logfc_res = adata.uns["logfc_res"]
    logfc_res.to_csv("scores/log2FC_LR_scores.tsv", sep="\t")
    # Make dotplot
    logfc_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="lr_logfc",
        size="lr_logfc",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=False,
        top_n=40,
        orderby="lr_logfc",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(12, 12),
        # uns_key to use, default is "liana_res"
        uns_key="logfc_res",
    )
    logfc_lr_dotplot.save(
        filename="figures/log2FC_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    logfc_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="means",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        top_n=40,
        orderby="lr_logfc",
        orderby_ascending=False,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="logfc_res",
        figure_size=(12, 12),
    )
    logfc_lr_tileplot.save(
        filename="figures/log2FC_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_cellchat(
    adata,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
    tileplot_width,
    tileplot_height,
):
    print(
        "\nPerforming LR analysis using CellChat method. Check liana.log for more information."
    )
    logging.info(cellchat.describe())
    cellchat(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        base=2,
        verbose=True,
        key_added="clcht_res",
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    clcht_res = adata.uns["clcht_res"]
    clcht_res.to_csv("scores/CellChat_LR_scores.tsv", sep="\t")
    # Make dotplot
    clcht_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="lr_probs",
        size="cellchat_pvals",
        # we inverse sign since we want small p-values to have large sizes
        inverse_size=True,
        top_n=40,
        orderby="cellchat_pvals",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(12, 12),
        # we filter the pvals column to <= 0.05
        filter_fun=lambda x: x["cellchat_pvals"] <= 0.05,
        # uns_key to use, default is "liana_res"
        uns_key="clcht_res",
    )
    clcht_lr_dotplot.save(
        filename="figures/CellChat_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )
    clcht_lr_tileplot = li.pl.tileplot(
        adata=adata,
        fill="trimean",
        label="props",
        label_fun=lambda x: f"{x:.2f}",
        filter_fun=lambda x: x["cellchat_pvals"] <= 0.05,
        top_n=40,
        orderby="cellchat_pvals",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        uns_key="clcht_res",
        figure_size=(12, 12),
    )
    clcht_lr_tileplot.save(
        filename="figures/CellChat_LR_TilePlot.pdf",
        width=tileplot_width,
        height=tileplot_height,
    )


def run_rankaggregate(
    adata,
    methods,
    celltype_column,
    resource_name,
    de_method,
    source_celltypes,
    target_celltypes,
    dotplot_width,
    dotplot_height,
):
    print(
        "\nPerforming LR analysis using Liana's implementation of aggregating all methods requested. Check liana.log for more information."
    )
    custom_rank_aggregate = li.mt.AggregateClass(li.mt.aggregate_meta, methods=methods)
    logging.info(custom_rank_aggregate.describe())
    custom_rank_aggregate(
        adata,
        groupby=celltype_column,
        expr_prop=0.1,
        resource_name=resource_name,
        verbose=True,
        n_perms=None,
        de_method=de_method,
    )
    # write the CellPhoneDB output to file
    liana_res = adata.uns["liana_res"]
    liana_res.to_csv("scores/CustomRankAggregate_LR_scores.tsv", sep="\t")
    # Make dotplot
    liana_lr_dotplot = li.pl.dotplot(
        adata=adata,
        colour="magnitude_rank",
        size="magnitude_rank",
        inverse_size=True,
        inverse_colour=True,
        top_n=40,
        orderby="magnitude_rank",
        orderby_ascending=True,
        source_labels=source_celltypes,
        target_labels=target_celltypes,
        figure_size=(8, 7),
        uns_key="liana_res",
    )
    liana_lr_dotplot.save(
        filename="figures/CustomRankAggregate_LR_DotPlot.pdf",
        width=dotplot_width,
        height=dotplot_height,
    )


def print_citations_to_helpfile():
    # ----- Print citations ------#
    citation1 = "Dimitrov D., Schäfer P.S.L, Farr E., Rodriguez Mier P., Lobentanzer S., Dugourd A., Tanevski J., Ramirez Flores R.O. and Saez-Rodriguez J. 2023 LIANA+: an all-in-one cell-cell communication framework. BioRxiv. https://www.biorxiv.org/content/10.1101/2023.08.19.553863v1"
    citation2 = "Dimitrov, D., Türei, D., Garrido-Rodriguez M., Burmedi P.L., Nagai, J.S., Boys, C., Flores, R.O.R., Kim, H., Szalai, B., Costa, I.G., Valdeolivas, A., Dugourd, A. and Saez-Rodriguez, J. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0"
    logging.info(
        f"\nUse the following citations for using Liana: \n{citation1}\n{citation2}\n"
    )


def provide_help(adata):
    print("\nCheck the liana_help.out file for more help running this script\n")
    with open("liana_help.out", "w") as helpfile:
        # ----- Some help copied from the liana-py.readthedocs.io ------#
        print("\nAvailable methods are:", file=helpfile)
        print(li.mt.show_methods(), file=helpfile)
        print("\nAvailable LR resources are:", file=helpfile)
        print(li.resource.show_resources(), file=helpfile)
        print("\nAvailable metadata in the AnnData object are:", file=helpfile)
        print(adata.obs.head(), file=helpfile)


def main():
    # ----- Get arguments from command line call ------#
    args = get_args()
    input_file = args.input_file
    provide_more_help = args.provide_more_help
    celltype_column = args.celltype_column
    source_celltypes = (
        None
        if ((args.source_celltypes is None) or (args.source_celltypes == "NULL"))
        else args.source_celltypes.split(sep=",")
    )
    target_celltypes = (
        None
        if ((args.target_celltypes is None) or (args.target_celltypes == "NULL"))
        else args.target_celltypes.split(sep=",")
    )
    resource_name = args.resource_name
    methodsToRun = args.methodsToRun.replace(" ", "").split(sep=",")
    de_method = args.de_method
    dotplot_width = args.dotplot_width
    dotplot_height = args.dotplot_height
    tileplot_width = args.tileplot_width
    tileplot_height = args.tileplot_height

    print_citations_to_helpfile()

    # ----- Load data ------#
    if (input_file in ["test_data", "NULL", ""]) or (input_file is None):
        print("\nNo input_file provided. Using pbmc68k_reduced testdata")
        print("And also setting some default values")
        adata = sc.datasets.pbmc68k_reduced()
        celltype_column = "bulk_labels"
        source_celltypes = ["CD34+", "CD56+ NK", "CD14+ Monocyte"]
        target_celltypes = ["CD34+", "CD56+ NK"]

    else:
        # load a h5ad file
        adata = sc.read_h5ad(input_file)

    # ----- Print help if that is what was asked for and quit ------#
    if provide_more_help is not None:
        provide_help(adata)
        quit()

    # Note: Liana uses .X count matrix by default. If this is not log1p normalized, make sure it is using the scanpy method:
    # scanpy.pp.log1p()
    # create umap with celltype_column annotations as a sanity check
    sc.pl.umap(adata, color=celltype_column, save=f"_{celltype_column}", show=False)
    print("\nFirst save UMAP plot")

    # Prepare folders for saving output
    os.makedirs(name="figures", exist_ok=True)
    os.makedirs(name="scores", exist_ok=True)
    # ----- The real analysis starts here ------#
    # Run the methods asked for and create some basic plots and export the output

    if "GeometricMean" in methodsToRun:
        run_geometric_mean(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    if "CellPhoneDB" in methodsToRun:
        run_cellphonedb(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    if "Connectome" in methodsToRun:
        run_connectome(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    if "log2FC" in methodsToRun:
        run_logfc(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    if "NATMI" in methodsToRun:
        run_natmi(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    if "SingleCellSignalR" in methodsToRun:
        run_singlecellsignalr(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    if "CellChat" in methodsToRun:
        run_cellchat(
            adata,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
            tileplot_width,
            tileplot_height,
        )

    # Convert names provided by the user to the liana methods
    methods_dict = {
        "CellPhoneDB": cellphonedb,
        "Connectome": connectome,
        "log2FC": logfc,
        "NATMI": natmi,
        "SingleCellSignalR": singlecellsignalr,
        "CellChat": cellchat,
        "GeometricMean": geometric_mean,
    }
    methods = [methods_dict[k] for k in methodsToRun if k in methods_dict.keys()]
    if len(methods) > 1:
        run_rankaggregate(
            adata,
            methods,
            celltype_column,
            resource_name,
            de_method,
            source_celltypes,
            target_celltypes,
            dotplot_width,
            dotplot_height,
        )


if __name__ == "__main__":
    main()
