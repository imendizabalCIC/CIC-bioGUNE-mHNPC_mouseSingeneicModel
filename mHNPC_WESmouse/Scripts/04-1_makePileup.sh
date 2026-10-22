#!/bin/bash

################################################################################
#| SAMTOOLS MPILEUP GENERATION - MOUSE WES
################################################################################
#| Last update: 22/10/2026
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#| This script generates a compressed pileup file from an aligned and
#| processed WES mouse BAM file using samtools mpileup.
#|
#| Main tasks:
#|   - Load the required Conda environment.
#|   - Extract read coverage and nucleotide information from the BAM file.
#|   - Restrict analysis to targeted exome regions defined in a BED file.
#|   - Compress and save the resulting pileup file.
################################################################################

#SBATCH --job-name=mpileup
#SBATCH -o logs/mpileup.out
#SBATCH -e logs/mpileup.err
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --partition=NORMAL
#SBATCH --mail-user=irondon@cicbiogune.es
#SBATCH --mail-type=FAIL

source /opt/ohpc/pub/apps/anaconda3/cic-env
conda activate WESenv

samtools mpileup -f ~/Sequences/Mus_musculus.GRCm38.dna.primary_assembly.fa -l ~/WES_data/mm10/Twist_Mouse_Exome_Chr_Modified_for_CODEX2.bed /vols/GPArkaitz_bigdata/irondon/AC-82_WESmouse_pipeline/02_Mapping/02_Mapped_Sequences/B844_mapped.rg.sorted.dup.bqsr.bam | gzip -c >myPileup.pileup.gz


git rm --cached mHNPC_WESmouse

mv mHNPC_WESmouse/.git ../mHNPC_WESmouse_git_backup

git add mHNPC_WESmouse/