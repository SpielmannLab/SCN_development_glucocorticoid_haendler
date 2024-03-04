# analysis_pipeline - readme

This pipeline starts from the initial seurat object created using the [count_to_seurat](../preprocessing/count_to_seurat/count_to_seurat_main.sh) script. The first step is to perform QC and clustering for a single-sample at a time using the [sc_per_sample](analysis_pipeline/1_sc_per_sample) pipline. Then merge the four samples using the [2_sc_multi_sample](analysis_pipeline/2_sc_multi_sample), perfom clustering and differential expression analysis. At the end, perform trajectory analysis using [3_sc_trajectory](analysis_pipeline/3_sc_trajectory). Cell-wise gene-set enrichment analysis was done using the scripts in [ssGSEA](ssGSEA/ssGSEA_sbatch.sh). Ligand-receptor interactions were analysed using the scripts in [cell_cell_signalling](cell_cell_signalling)  In all cases, read the instructions at the *_main.sh, *_sbatch.sh, or the associated readme.md files.

## Requirements

This script relies on the existing installation of nextflow (version 22.04.1) via [modules](https://modules.readthedocs.io/en/latest/) as well as the conda environments *scVelocity* and *scTrajectory* as mentioned in the main installation instructions.
