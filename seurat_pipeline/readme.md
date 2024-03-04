# seurat_pipeline - readme

This pipeline starts from the initial seurat object created using the [count_to_seurat](../preprocessing/count_to_seurat/count_to_seurat_main.sh) script. The first step is to perform QC and clustering for a single-sample at a time using the [sc_per_sample](seurat_pipeline/1_sc_per_sample) pipline. Then merge the four samples using the [2_sc_multi_sample](seurat_pipeline/2_sc_multi_sample), perfom clustering and differential expression analysis. At the end, perform trajectory analysis using [3_sc_trajectory](seurat_pipeline/3_sc_trajectory). In all cases, read the instructions at the *_main.sh script.

## Requirements

This script relies on the existing installation of nextflow (version 22.04.1) via [modules](https://modules.readthedocs.io/en/latest/) as well as the conda environments *scVelocity* and *scTrajectory* as mentioned in the main installation instructions.
