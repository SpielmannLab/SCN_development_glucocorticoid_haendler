#! /bin/bash

### Submit this Script with: sbatch <script.sh> ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition shortterm/debug/longterm:
#SBATCH --partition=shortterm
#  Use so many node:
#SBATCH --nodes=1
#  Request so many cores (hard constraint):
#SBATCH -c 3
#  Request so much of memory (hard constraint):
#SBATCH --mem=350GB
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

# Load your necessary modules:
module load nextflow/v22.04.1

# Move to SCRATCH were everything will be done
cp *.nf $SCRATCH/ && \
    cp *.yaml $SCRATCH/ && \
    cp -r bin $SCRATCH/ && \
    cd $SCRATCH && \
    chmod +x bin/*

tree

# Submit the Nextflow Script:
nextflow run ssGSEA.nf -params-file ssGSEA_params.yaml --id $SLURM_JOB_ID -resume
