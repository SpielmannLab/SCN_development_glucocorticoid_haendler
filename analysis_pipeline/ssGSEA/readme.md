# Installation of conda environments for ssGSEA pipeline
*Note* - the installations can be sped up by using "mamba" instead of "conda", but it needs to be installed first.

- For **ssGSEA** analysis, a separate conda environment called "ssGSEA" is necessary. Do either of the following to set it up:

    - Option 1 (takes long about 2-3 days, but portable to other HPCs)

            # ssh into the omics cluster
            srun --partition=debug --mem=50GB -c 1 --pty bash
            conda create --name ssGSEA "conda-forge::r-seurat=4.1.1" "conda-forge::r-ggplot2>3.0" "conda-forge::r-tidyverse=1.3.2" "bioconda::bioconductor-gseabase" "bioconda::bioconductor-gsva" "conda-forge::r-msigdbr" "bioconda::bioconductor-biocparallel" "conda-forge::r-matrix" "conda-forge::r-reshape2" "conda-forge::r-stringr" "bioconda::bioconductor-matrixgenerics" "conda-forge::r-broom" "conda-forge::r-biocmanager" "conda-forge::r-patchwork" "bioconda::bioconductor-limma" "bioconda::bioconductor-biocversion"
            
            # Install additional packages
            conda activate ssGSEA
            conda install -c conda-forge r-docopt
            R
            BiocManager::install("escape", update=FALSE)
            install.packages("docopt")
            install.packages("pheatmap")
        
    - Option 2 (much faster, but limited to Uni-Luebeck OMICS cluster)

            # ssh into the omics cluster
            srun --partition=debug --mem=50GB -c 1 --pty bash
            conda env create --file installation/ssGSEA.yml 
            
            # Install additional packages
            conda activate ssGSEA
            R
            BiocManager::install("escape", update=FALSE)
            install.packages("docopt")
            install.packages("pheatmap")

- Statistics and Violin plots in sGSEA requires another conda environment, called "ggplot_essentials"
