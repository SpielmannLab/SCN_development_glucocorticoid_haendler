# installation - readme

All the files with *.yml extensions contain conda environment recipes. These environments can be installed using [mamba](https://github.com/conda-forge/miniforge) or conda using the following command (example):

    mamba env create -f liana_dev.yml

Note, for the conda environment scVelocity, after creation of the conda environment using [scVelocity.yml](installation/scVelocity.yml), it is necessary to run the [install_packages.R](installation/install_packages_in_scVelocity.R) from within the scVelocity conda environment to install additional R-packages.

Note, for the conda environment ssGSEA, please refer to [the specfic readme file](../analysis_pipeline/ssGSEA/readme.md).

## Installation time

The installation of all these environments should take less than 10 minutes using mamba. However, the installation of ssGSEA.yml can take much longer (see the corresponding readme file). The installation of R packages within the scVelocity and ssGSEA environments (see above) can take another 10 minutes.
