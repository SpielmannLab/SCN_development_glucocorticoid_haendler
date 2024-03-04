# Cell cell signalling

The workflow allows two specific cell-cell-signalling tools to be run

1. NICHES - in differential mode only
2. LIANA - in single-condition as well as differential mode

## Installation requirements

This workflow is self contained. However it requires that you have mamba installed. This is a faster implementation of conda. If you do not have it set up, you can use just do `conda install mamba` from your base conda environment.

## Important considerations

1. Differential NICHES (Cell2Cell) requires that the number of cells in the two conditions be roughly the same.
2. Only Cell2Cell mode of analysis has been tested with NICHES. Other possible modes System2Cell and Cell2System.
3. Differential LIANA requires that each condition has multiple samples. This is because py-DESeq2 is used to calculate the differentially expressed genes.

## Expected output files and their Interpretation

1. NICHES

   The results are that from differential NICHES analysis between two conditions. In NICHES the LR interactions that may be occuring in this dataset are quantified on a cell-cell pair basis (i.e., not cluster-cluster basis). The pipeline generates the three main outputs as well as a raw data file:

   - A 3-panel UMAP plot, the dots represent cell-cell pairs and the dimensions are strength of LR signalling (instead of gene expression in a normal scRNA-seq UMAP plot). While the UMAP is the same across the three panels, their coloring is different:

     - Top panel: Colored based on the cell types of the Sending-Receiving cells
     - Middle panel: Colored based on the cell types of the Sending cells
     - Bottom panel: Colored based on the cell types of the Receiving cells

   - A Volcano plot, where the signalling strength between cells for each LR gene-pair is compared between the two conditions. Therefore, LR pairs are the dots and colors are celltype-celltype pairs. That is points to the extreme left and right in the volcano plot that pass the p-value significance threshold are likely to be ligand-receptor pairs that differentially signalling between the two conditions in those sending-receiving cell type pair. There are again three panels:

     - Top-left panel: Colored based on the cell types of the Sending cells
     - Top-right panel: Colored based on the cell types of the Receiving cells
     - Bottom panel: Colored based on the cell types of the Sending-Receiving cells

   - A TSV file, containing the differential cell-cell-signalling information summarized in the Volcano plot.
   - An RDS file (containing the cell-cell signalling Seurat object) containing the processed data.

2. LIANA

   LIANA is a meta analysis tool that performs many different cell-cell-signalling analyses under the hood per sample condition. Unlike NICHES, all the quantifications occur at the cell-cluster (i.e., cell-type) level. The following two figure files are generated per algorithm.

   - A tile plot: This plot shows the level of ligand and receptor expression for the top 40 LR pairs. This color bar labeled in this plot shows the expression level of the ligand-gene in the source cell type (i.e., sender cell type) and the expression level of the receptor-gene in the target cell type (i.e., the receiver cell type). The top 40 LR pairs are selected based on how specific that signalling mechanism is across all the possible sender-receiver celltype combinations. The only exception to this ordering is for the tileplot from SingleCellSignalR algorithm, where it is ordered based on the magnitude of LR interaction.

   - A dot plot of LR magnitude, specificity and expression: This plot shows the the magnitude and specificity of LR interactions between different cell types. In most algorithms, the magnitude is used to color the dots (see scale bar - higher the number higher the interaction), and the size of thd dot represents the specificity (bigger the dot, more specific is this interaction between those two cell types - note: -log10 has been applied to p-values). Some exceptions include the "CustomRankAggregate" and "SingleCellSignalR" methods, which do not provide any specificity values and the "log2FC" method, which does not provide any magnitude values.
   - Scores for downstream analysis
   -

3. Differential LIANA (diff_liana)
   This workflow, which analysis the differentially regulated ligand-receptor pairs generates the following outputs. This workflow is not optimised yet. But the idea is that only those LR pairs are plotted here, for which the Ligand, the Receptor-related genes, and the LR-combined are significantly different between the two conditions.

   - Tile plots of the differentially regulated ligands and receptors (see above for a description of tile plot)
   - QC of pseudobulk

## Regulon specific analysis

To restrict the tile plot of differentially expressed ligands and receptors to only those genes that are in the Nr3c1-regulon, the following scripts were used:

    '''R
    # Collection of regulons for Nr3c1 from collecTRI
    set <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE, ..)
    set_interest <- set %>% dplyr::filter(source == "Nr3c1")
    write.table(set_interest, file = "collectTRI_Nr3c1_regulons.tsv", sep = "\t", col.names = NA)
    set_interest_activators <- set_interest %>% dplyr::filter(mor > 0)
    write.table(set_interest_activators, file = "collectTRI_Nr3c1_regulons_activating.tsv", sep = "\t", col.names = NA)
    set_interest_repressors <- set_interest %>% dplyr::filter(mor < 0)
    write.table(set_interest_repressors, file = "collectTRI_Nr3c1_regulons_repressing.tsv", sep = "\t", col.names = NA)
    '''

Then the differentially expressed genes obtained from the output *diff_lian/xx_scn_scores/* directory was filtered against the above gene list using the script [filter_diff_liana_lr_for_geneset.R](utils/filter_diff_liana_lr_for_geneset.R) and plotted using the script [make_tileplots_for_filtered_diff_liana_LR.py](utils/make_tileplots_for_filtered_diff_liana_LR.py) 
