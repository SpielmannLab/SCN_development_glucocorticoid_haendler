#! /bin/bash
# Run this script to get cell x gene matrix from fastq files
# STEP1. Make sure you have run the mkfastq_main.sh script
# STEP2. Provide all the details in the INPUT parameters section
# STEP3. Close this file and execute it using ./count_main.sh


# INPUT Parameters
# This is where the all the fastq files are located "/*/*/*,/*/*/*,/*/*"
path_fastqdata="/data/SCN/fastq"
# Path of the reference transcriptome /*/*/*
path_refdata="/data/refdata-gex-mm10-2020-A"
# This is the place where output will be stored. The directory should exist. /*/*/*
path_out="/data/SCN/counts"
# Samples name. Eg. For fastq files of the format "MySample_S1_L001_R1_001.fastq.gz" use "MySample"
# Provide multiple samples as ( "MySample1" "MySample2" "MySample3" )
# IMPORTANT: The fastq files corresponding to these samples must be present in all the path_fastqdata provided. Else run the samples separately.
samples=( "GD17" "PND02" "PND10" "PND30")
# Whether or not to include reads mapped to introns in the count matrix
intron=true # has to be lower case true or false 

# DO NOT CHANGE THINGS HERE
for sample in ${samples[@]}
do
	if [[ -d ${path_out}/${sample} ]]; then
		echo "Error: The directory ${path_out}/${sample} already exists"
		exit 1
	else
		sbatch -c 16 --job-name=count_${sample} count_core.sh $sample $path_refdata $path_fastqdata $path_out $intron
		echo "Message: Job count_${sample} sent to count pipline"
		echo "Look for slurm-count_${sample}-xxxxxx.out to follow progress"
	fi
done
