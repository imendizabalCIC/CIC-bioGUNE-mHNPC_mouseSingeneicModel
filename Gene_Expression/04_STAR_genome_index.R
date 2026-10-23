#################################################################################
#                             STAR GENOME INDEX
#################################################################################

# Summary
#---------

# The reads will be mapped using STAR which is a genome reference based aligner. 
# To run STAR two steps must be followed:
#  - Step 1. Build STAR genome index (00_STAR_genome_index.sh)
#  - Step 2. Mapping reads with STAR is performed in this script
#
# In this script, the first step is performed which is building the genome index 
# based on the insert length. This can only be done for human and mouse. 
# 
# This step is not mandatory every time you run STAR. The multiple genome indexes 
# have being previously done which means you just have to load the file. 





#################################################################################
#                           STAR VERSION
#################################################################################


# STAR version 2.7.10a

## Manuals used 

# STAR manual 2.7.10a
# October 20, 2022
# Link: https://github.com/alexdobin/STAR

# STAR manual 2.7.0a
# January 23, 2019
# Link: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf


#################################################################################
#                                 PIPELINE
#################################################################################

### Access R ###
R
# Load R libraries 
.libPaths(".libPaths")


#-------------------------------------------------------------------------------------------------------------------------------------------------------
### General Project ###

# Pathway to the folders and files
# Select one option depending if you are running the script in Rocky or local
path <- "path"

# Specie
# specie <- "Mouse"
specie <- "Human"

# Read insert
# Remove one nucleotide
# read <- 100
read <- 150

### Process information ###
partition <- "FAST"
# time <- c("00:45:00") # Mouse genome
time <- c("02:00:00") # Human genome
memory <- c("40")
cpu <- 8
ram <- 30000000000 
#-------------------------------------------------------------------------------------------------------------------------------------------------------

# Input and output directories  
indexfolder <- "Indexes"
folder_name <- paste(specie, read+1, sep = "_")
dir.create(file.path(indexfolder, folder_name))
dir_outfiles <- paste(indexfolder,"/", folder_name, sep = "")

# Set directory
setwd(indexfolder)

# Select specie annotation files
if(specie == "Human"){
  # Human
  fasta <- "Homo_sapiens.GRCh38.dna.primary_assembly_110.fa" 
  gtf <- "Homo_sapiens.GRCh38.110.gtf"
}else{
  # Mouse
  fasta <- "Mus_musculus.GRCm39.dna.primary_assembly_110.fa" 
  gtf <- "Mus_musculus.GRCm39.110.gtf"
}


###############################################################
#                           STAR
###############################################################

### Generate PBS ###

job_name <-  paste("Genome_index", specie, read, format(Sys.time(),"%Y%m%d"), sep = "_")

# STAR path
star_path <- "source environment"

# Command
command <- paste("STAR --runThreadN", cpu,    
                 "--runMode genomeGenerate",
                 "--genomeDir", dir_outfiles,                     
                 "--genomeFastaFiles", fasta,              
                 "--sjdbGTFfile", gtf,
                 "--sjdbOverhang", read, 
                 sep=" ") 

# SBATCH File
filename <- paste(job_name,".sh",sep='');
cat(
  c("#!/bin/sh"),
  c("#SBATCH  --export=ALL"),
  paste("#SBATCH --job-name=",job_name,sep=''),
  paste("#SBATCH --partition=",partition,sep=''),
  paste("#SBATCH --cpus-per-task=",cpu,sep=''),
  paste("#SBATCH --time=",time,sep=''),
  paste("#SBATCH --mem=",memory,"GB",sep=''),
  paste("#SBATCH -o ",job_name,".out",sep=''),
  paste("#SBATCH -e ",job_name,".err",sep=''),
  "\n\n",
  c(star_path), 
  "\n",
  c(paste("cd ", indexfolder)), 
  "\n",
  c(command),
  file=filename,sep = "\n",append=F)

# Run PBS
system(paste("sbatch",filename,sep=' '));
Sys.sleep(3)

q()
