# Readme

These scripts were used to pre-process the raw seqencing data, and analyse the resultant count matrix, as well as to prepare figures in the manuscript: **[Development of the central circadian clock sensitivity to glucocorticoids](bioRxiv link)

## Main steps
Here is a summary of all the analysis steps.

The scripts for each step are contained within the respective directories and an associated readme file, which contain individual installation and usage requirements, whenever necesary. Down below are some global installation instructions and system requirements.

1. [Preprocessing](preprocessing)
10X Cellranger-based scripts for obtaining single-cell gene-level count matrixes from fastq files.

2. [seurat_pipeline]](seurat_pipeline)
Contains generic NextFlow-based pipeline for basic sn-RNAseq data analysis that includes QC reporting, QC filtering including doublet detection, clustering, and marker-gene identification, merging/integration of multiple samples, and trajectory analysis. These are implemented using Seurat (v4) and monocle3, in addition to other R-packages.

3. **custom_analysis** contains scripts developed to perform custom analysis and custom figures for this work.

## Installation instructions
Most of these scripts were written to work in a HPC running Debian GNU/Linux 11 (bullseye). Most scripts are wrapped within an additional [Slurm](https://slurm.schedmd.com) job-submission script (sbatch), which can be bypassed if running on a local machine. Moreover, all the paths needs to be adjusted as required.

All the packages were installed within conda environments. Instructions/commands to generate these conda environments are located in installation/*.yaml

Note, for the conda environment scVelocity, after creation of the conda environment using [scVelocity.yml](installation/scVelocity.yml), it is necessary to run the [install_packages.R](instllation/install_packages.R) from within the scVelocity conda environment to install additional R-packages.
