#! /bin/bash
# Run this script to perform ssGSEA analysis
# STEP1. Provide all the details in the file "ssGSEA_params.yaml"
# STEP2. Close this file and execute it using ./ssGSEA_main.sh. You may have to run the command <chmod +x ssGSEA_main.sh> before hand.

sbatch --job-name=ssGSEA ssGSEA_core.sh
