# installation - readme

All the files with *.yml extensions contain conda environment recipes. These environments can be installed using [mamba](https://github.com/conda-forge/miniforge) or conda using the following command (example):

    mamba env create -f liana_dev.yml

Note, for the conda environment scVelocity, after creation of the conda environment using [scVelocity.yml](installation/scVelocity.yml), it is necessary to run the [install_packages.R](installation/install_packages_in_scVelocity.R) from within the scVelocity conda environment to install additional R-packages.
