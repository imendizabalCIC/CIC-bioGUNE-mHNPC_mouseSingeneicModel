#!/bin/bash

#################################################################################
#| QUALITY CONTROL ON THE FASTQ FILE WITH FASTQC SOFTWARE FOR RNASEQ ANALYSIS ##
#|    READ 2
#################################################################################
#| Date: 22/04/2024
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This SLURM job runs FastQC on all Read 2 FASTQ files from the AC-45_RNAseq-FFPE 
#| project. FastQC generates quality control reports (HTML + zipped folders) for 
#| each sample, assessing per-base quality, GC content, adapter contamination, and 
#| other metrics. 
#|
#| Workflow:
#|   - Input: Trimmed Read 2 FASTQ files (*.fastq.gz) from $FASTQs_dir
#|   - Tool:  FastQC (run with 2 CPUs)
#|   - Output: Quality reports stored in $FASTQs_out
#|
#| Notes:
#|   - Only Read 2 files are analyzed in this script
#|   - Folder $FASTQs_out is created if it does not already exist.
##################################################################################


#SBATCH --job-name="FASTQC_READ_2"
#SBATCH -o FASTQC_READ2.out
#SBATCH -e FASTQC_READ2.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=3GB
#SBATCH --partition=FAST
#SBATCH --time=12:00:00

source /opt/ohpc/pub/apps/anaconda3/cic-env 
conda activate WESenv

#| Defining the working directory
projectDir="/AC-82_WESmouse/"
numProcessors=2
FASTQs_out="/sharedir/AC-82_WESmouse/FASTQCs"
FASTQs_dir="/sharedir/AC-82_WESmouse/FASTQs_trimmed"
cd $projectDir

#| Create folder to store the output values. If the folder exist it does not modify the folder
mkdir $FASTQs_out 

#| Quality control of the sequences from FASTqs files (all the files in the folder)
fastqc $FASTQs_dir/*_2.fastq.gz -o $FASTQs_out/ -t $numProcessors