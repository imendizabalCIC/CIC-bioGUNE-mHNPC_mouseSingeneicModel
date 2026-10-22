#!/bin/bash

################################################################################
#| COPY NUMBER ALTERATION ANALYSIS - CONTROL-FREEC - MOUSE WES
################################################################################
#| Last update: 22/10/2026
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script submits a SLURM job to detect copy number alterations (CNAs)
#| in WES mouse samples using Control-FREEC.
#|
#| Main tasks:
#|   - Allocate computational resources through SLURM.
#|   - Activate the required Conda environment.
#|   - Run Control-FREEC using a sample-specific configuration file.
#|   - Generate CNA results for downstream analysis.
################################################################################

#SBATCH --job-name=J993
#SBATCH -o logs/J993_Control-FREEC.out
#SBATCH -e logs/J993_Control-FREEC.err
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --partition=NORMAL
#SBATCH --mail-user=irondon@cicbiogune.es
#SBATCH --mail-type=FAIL

source /opt/ohpc/pub/apps/anaconda3/cic-env
conda activate /vols/GPArkaitz_bigdata/DATA_shared/NewCluster_Software/conda_envs/WESenv

freec -conf /vols/GPArkaitz_bigdata/irondon/AC-82_WESmouse_pipeline/05_CNAs_CONTROL_FREEC_MP/J993_config_exome.txt
