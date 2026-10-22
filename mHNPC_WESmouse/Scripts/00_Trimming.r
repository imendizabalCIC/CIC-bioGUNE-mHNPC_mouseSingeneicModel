#/opt/R/R-4.1.2/bin/R

###############################################################################
#| THIS CODES CREATES A SCRIPT FOR EACH FASTQ FILE AND SUMIT JOBS 
#| TO PERFORM TRIMMING ON THE FASTQ FILES OF THE SECOND READ      
###############################################################################
#| Date: 22/04/2024
#| Author: Ivana Rondon Lorefice
#|
#| This script creates and execute jobs for trimming corresponding DNA sequences 
#|
#| The raw data analyzed is contained in the folder /sharedir/AC-82_WESmouse/FASTQs/
#|
#| Library specification:
#| Libraries were sequenced on the NovaSeq 6000 platform, originating paired-end reads of 150 base pairs in length, 
#| with an average of 113 million reads per sample. Adaptor sequences specified by the Twist Universal Adapter System - TruSeq Compatible
#|
#| In summary: This scripts creates a .sh file and then it is launched to the job 
#| executor in indar to:
#|   1) Cut the Twist Universal Adapter System - TruSeq Compatible illumina adaptors from the fastq files
#|   2) Perform quality trimming with a cutoff of 10
#|   3) Minimum length size of 20bp. This is chosen based on the distribution of 
#| length for our insert.
#|
#| From the cutadapt documentation (paired-end sequence)
#|     cutadapt [options] -o output 1 -p output2 input1 input2
#| [options]: 
#|     -a: Cut adapters of the 3' of Read 1
#|     -A: Cut adapters of the 3' of Read 2
#|     -m: Defines a threshold of the minimum lenght sequence to keep
#|     -q: Defines a threshold of the quality for each sequence
#|     -u: Cut a specific number of nucleotides of Read 1
#|     -U: Cut a specific number of nucleotides of Read 2
###############################################################################


################################################################################
#|  DATA LIST  
################################################################################
# List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = "/sharedir/AC-82_WESmouse/FASTQs/", pattern = "*_2.fastq.gz")
samples <- gsub("_2.fastq.gz", "", samples)
################################################################################


################################################################################
#|  DEFINED PARAMETERS
################################################################################
# Defining variables
cpu <- "2"
memory <- "1G"
time <- "5:00:00"
partition <- "FAST"
################################################################################


################################################################################
#|  JOBS LOOP  
################################################################################
#| This creates a loop over all the samples and create a .sh file which is then submitted as a job

for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste("TRIMMING_", samples[s], ".sh", sep = "") 

  # Forward strand  (Read 1)
  input1 <- paste("/sharedir/AC-82_WESmouse/FASTQs/", samples[s], "_1.fastq.gz", sep = "")
  input2 <- paste("/sharedir/AC-82_WESmouse/FASTQs/", samples[s], "_2.fastq.gz", sep = "")

  # Output directories
  output1 <- paste("/sharedir/AC-82_WESmouse/FASTQs_trimmed/", samples[s], "_trimmed_1.fastq.gz", sep = "") 
  output2 <- paste("/sharedir/AC-82_WESmouse/FASTQs_trimmed/", samples[s], "_trimmed_2.fastq.gz", sep = "")
  
  cat(
    "#!/bin/sh",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o ", samples[s], "_TRIM.out", sep = ""),
    paste("#SBATCH -e ", samples[s], "_TRIM.err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    "#SBATCH --mail-user=irondon@cicbiogune.es",
    "#SBATCH --mail-type=FAIL",
    "\n\n",
    "source /opt/ohpc/pub/apps/anaconda3/cic-env", 
    "conda activate cutadaptenv",
    "\n",
    paste("cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 0 -q 10 -m 20 --pair-filter=any -o", output1, "-p", output2, input1, input2, sep = " "),
    file = filename, sep = "\n", append = FALSE)
  system(paste("sbatch", filename, sep = " "))
}
