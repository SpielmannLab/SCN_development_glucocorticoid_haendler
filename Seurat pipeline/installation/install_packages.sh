#! /bin/bash
### Submit this Script with: sbatch script.sh ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=debug
#  Use one node:
#SBATCH --nodes=1
#  Request 10 cores (hard constraint):
#SBATCH -c 2
#  Request 64GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Find your job easier with a name:
#SBATCH --job-name=InstallingPackages

# Load your necessary modules
PATH=$WORK/.omics/anaconda3/bin:$PATH
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh
conda activate scVelocity

Rscript /data/humangen_mouse/scpipeline/installation/install_packages.R 
