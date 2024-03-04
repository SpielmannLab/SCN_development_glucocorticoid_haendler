# pre_processing - readme

The pre_processing is comprised of the following two steps:

## Creating BAM file and counting gene expression using Cellranger

Run this script to get cell x gene count matrix from fastq files. Follow the instructions in [count_main.sh](cellranger/count_main.sh) script. **Input**: Provide the path to the fastq files, reference transcriptome in the [count_main.sh](cellranger/count_main.sh).

### Requirements for creating gene expression matrix

This script relies on the existing installation of 10x cellranger version 5.0.1 via [modules](https://modules.readthedocs.io/en/latest/) and submits a sbatch job for every sample provided.

## Creating Seurat object (*.rds* file) from count matrix

This is a nextflow based script to convert count-matrix (output of cellranger - count pipeline) into Seurat object, using the *Read10X* function.  **Optionally**, it also runs the velocyto pipeline to create a loom file. The spliced and unspliced count matrices from the loom file is then merged with the Seurat object. Follow the instructions in [count_to_seurat_main.sh](count_to_seurat/count_to_seurat_main.sh). **Input**: Requires several files from cellranger.

### Requirements for creating Seurat object

This script relies on the existing installation of nextflow (version 22.04.1) via [modules](https://modules.readthedocs.io/en/latest/) as well as the conda environment *scVelocity* as mentioned in the main installation instructions.
