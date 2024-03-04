packages <- c('ggplot2','cowplot','future','dplyr','ggwordcloud','limma','pcaMethods')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos="https://ftp.gwdg.de/pub/misc/cran/")
    BiocManager::install(packages)

install.packages('httr', repos="https://ftp.gwdg.de/pub/misc/cran/")
install.packages('plotly', repos="https://ftp.gwdg.de/pub/misc/cran/")
library(reticulate)
conda_install("r-reticulate", "pandas")
install.packages('remotes', repos="https://ftp.gwdg.de/pub/misc/cran/")
# Force installation of Seurat4
# install.packages('Seurat', repos="https://ftp.gwdg.de/pub/misc/cran/")
remotes::install_version("SeuratObject",
  "4.1.4",
  repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

install.packages('docopt', repos="https://ftp.gwdg.de/pub/misc/cran/")
install.packages('ggwordcloud', repos="https://ftp.gwdg.de/pub/misc/cran/")
install.packages('harmony', repos="https://ftp.gwdg.de/pub/misc/cran/")
install.packages('R.utils', repos="https://ftp.gwdg.de/pub/misc/cran/")
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github("mojaveazure/seurat-disk")
