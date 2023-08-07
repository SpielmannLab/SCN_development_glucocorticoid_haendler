#! /bin/bash    
PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

# Load your necessary modules:
module load nextflow/v22.04.1

# Move to SCRATCH were everything will be done
cp -Ru bin $SCRATCH/
cp -u *.yaml $SCRATCH/
cp -u *.nf $SCRATCH/

cd $SCRATCH
chmod +x bin/*

# Use conda environments
NXF_CONDA_ENABLED=true

# Submit the Nextflow Script:
nextflow run TOME.nf -params-file TOME.yaml --id ${SCRATCH/"/scratch/"/}