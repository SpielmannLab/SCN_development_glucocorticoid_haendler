# analysis_pipeline - readme

This pipeline starts from the initial seurat object created using the [count_to_seurat](../pre_processing/count_to_seurat/count_to_seurat_main.sh) script.

1. The first step is to perform QC and clustering for a single-sample at a time using the [sc_per_sample](analysis_pipeline/sc_per_sample) pipeline.

    - Expected time for the demo data: 20 minutes
    - Expected output: Multiple folders will be created at the end of the pipeline. These include QC plots, filtering statistics, clustering and marker genes along with seurat objects as \*.rds files. Use the dim_reduc/\*.rds file to continue.

2. Then merge the four samples using the [sc_multi_sample](analysis_pipeline/sc_multi_sample), perfom clustering and differential expression analysis.
    - Expected time for the demo data: 20 minutes
    - Expected output: Multiple folders will be created at the end of the pipeline. These include QC plots, filtering statistics, clustering and marker genes along with seurat objects as \*.rds files. Use cluster/\*.rds for trajectory analysis.

3. Perform trajectory analysis using [sc_trajectory](analysis_pipeline/sc_trajectory).
    - Expected time for the demo data: 5 minutes
    - Expected output: Three output directories will be created (pseudotime, pseudotime_in_seurat, and trajectory). The figures from pseudotime_in_seurat were used in the manuscript.

4. Cell-wise gene-set enrichment analysis was done using the scripts in [ssGSEA](ssGSEA/ssGSEA_sbatch.sh).
    - Expected time for the demo data: 30 minutes
    - Expected output: Multiple output directories will be created based on the parameters in the yaml file. These can include the raw scores in the gsea directory, in addition to the heatmap files in the heatmap directory.

5. Ligand-receptor interactions were analysed using the scripts in [cell_cell_signalling](cell_cell_signalling)  In all cases, read the instructions at the \*_main.sh, \*_sbatch.sh, or the associated readme.md files.
    - Expected time for the demo data: 10 minutes
    - Expected output: Multiple output directories are created based on the parameters in the yaml file. These can include NICHES, scanpy_objects, diff_liana, etc. The diff_liana/\*_scn_scores/Dysregulated_LRs_unfiltered.tsv contains the Ligand-receptor scores that were used to create the figure in the manuscript. To customise which LR-pairs are plotted, refer to the [readme file within the cell_cell_signalling](analysis_pipeline/cell_cell_signalling/README.md)

## Requirements

This script relies on the existing installation of nextflow (version 22.04.1) via [modules](https://modules.readthedocs.io/en/latest/) as well as several conda environments as mentioned in the main installation instructions.
