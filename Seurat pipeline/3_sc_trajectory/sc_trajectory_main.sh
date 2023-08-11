#! /bin/bash
# Run this script to perform trajectory analysis
# STEP1. Make sure you have an *.rds file containing a seurat object
# STEP2. Provide all the details in the file "sc_trajectory_params.yaml"
# STEP3. Close this file and execute it using ./sc_trajectory_main.sh

sbatch sc_trajectory_core.sh 
