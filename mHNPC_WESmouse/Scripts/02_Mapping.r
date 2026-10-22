#/opt/R/R-4.1.2/bin/R

################################################################################
#| MAPPING SEQUENCES READS
################################################################################
#| Date 09/04/2025
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#|
#| Generates and submits a per-sample SLURM job for NGS/WES preprocessing.
#| For each paired FASTQ, the job:
#|   - Aligns reads with BWA-MEM
#|   - Adds read groups (Picard)
#|   - Sorts & indexes BAM (samtools)
#|   - Marks/removes duplicates (Picard)
#|   - Performs BQSR (GATK BaseRecalibrator + ApplyBQSR)
#|   - Runs basic QC (samtools flagstat)
#| Uses environment variables (.env) to locate references, inputs, and outputs,
#| and to set SLURM resources (CPU, memory, time, partition, email).
#|
#| Workflow: 
#|
#| 1) Config & Discovery
#|    - Load .env
#|    - Set library paths
#|    - Find *_trimmed_1.fastq.gz and derive sample IDs
#|
#| 2) Per-Sample Job Script (02_Mapping_<sample>.sh)
#|    - Conda: source activate <WES env>
#|    - ALIGN:   bwa mem -t <cpu> <ref.fa> R1 R2 | samtools view -@ <cpu> -Sb > <sample>_mapped.bam
#|    - RG:      picard AddOrReplaceReadGroups I=<mapped.bam> O=<rg.bam> RGID=<id> RGLB=Library1 RGPL=illumina RGPU=Unit1 RGSM=<sample>
#|    - SORT:    samtools sort -@ <cpu> <rg.bam> -o <rg.sorted.bam> && samtools index <rg.sorted.bam>
#|    - DUPS:    picard MarkDuplicates INPUT=<rg.sorted.bam> OUTPUT=<rg.sorted.dup.bam> \
#|               METRICS_FILE=<metrics.txt> CREATE_INDEX=true REMOVE_DUPLICATES=true
#|    - BQSR-1:  gatk BaseRecalibrator -I <rg.sorted.dup.bam> -R <ref.fa> -L <intervals.list> \
#|               --known-sites <known.vcf[.gz]> -O <recal.table>
#|    - BQSR-2:  gatk ApplyBQSR -I <rg.sorted.dup.bam> -R <ref.fa> -L <intervals.list> \
#|               --bqsr-recal-file <recal.table> -O <rg.sorted.dup.bqsr.bam>
#|    - INDEX:   samtools index <rg.sorted.dup.bqsr.bam>
#|    - QC:      mkdir -p <mapping_qc_dir>/qualimap_<sample>
#|               samtools flagstat <rg.sorted.dup.bqsr.bam> > <mapping_qc_dir>/<sample>_samtools_qc.txt
#|               # (optional) qualimap bamqc ... ; coverage depth (bedtools) disabled by default
#|
#| 3) Submission
#|    - Write job script with SBATCH directives:
#|        --job-name=<sample>
#|        --cpus-per-task=<cpu>, --mem=<mem>, --time=<hh:mm:ss>, --partition=<part>
#|        -o logs/02_Mapping_<sample>.out, -e logs/02_Mapping_<sample>.err
#|        --mail-user=<email> --mail-type=FAIL
#|    - Submit via: sbatch 02_Mapping_<sample>.sh
#|
#| Outputs:
#|
#| <MappedDir>/<sample>_mapped.bam
#| <MappedDir>/<sample>_mapped.rg.bam
#| <MappedDir>/<sample>_mapped.rg.sorted.bam (+ .bai)
#| <MappedDir>/<sample>_mapped.rg.sorted.dup.bam (+ .bai)
#| <MappedDir>/<sample>_mapped.rg.sorted.dup.bqsr.bam (+ .bai)
#| <RecallDir>/recal_data_<sample>.table
#| <MarkedDupDir>/marked_dup_metrics_<sample>.txt
#| <MappingQCDir>/<sample>_samtools_qc.txt
#| logs/02_Mapping_<sample>.out, logs/02_Mapping_<sample>.err
#|
#| Comments:
#|
#| - Reference FASTA must be indexed for BWA (bwa index) and have accompanying .fai/.dict as needed.
#| - Known variant sites (SNP+indels) should be bgzipped & tabix-indexed for GATK.
#| - Duplicate removal is enabled (REMOVE_DUPLICATES=true); adjust per downstream needs.
#| - Qualimap and depth-by-bedtools are scaffolded but commented out for speed.
#| - Designed for mouse data but organism-agnostic with appropriate references.
#| - All paths/resources pulled from .env: BASE_DIR, DATA_SHARED, FASTQs_DIR, GENOME, SNPS_INDELS_DATA,
#|   INTERVAL_LIST, BED_FILE, directories for outputs, and SLURM resources.
################################################################################


################################################################################
#|  DATA CONFIGURATION  
################################################################################

# Configuration env file
conf_env <- "../00_conf_env.env"

# Loading necessary libraries and data
.libPaths("/R_libraries")
library("dotenv")
load_dot_env(file = conf_env)  # Carga las variables del archivo .env

# Setting based list directories
Base_dir <- Sys.getenv("BASE_DIR")
Data_shared <- Sys.getenv("DATA_SHARED")

# FASTQs directories
FATSQs_directory <- paste0(Data_shared, Sys.getenv("FASTQs_DIR"))

# List with the names of the files we will process, only keeping the sample name.
samples <- list.files(path = FATSQs_directory, pattern = "*_trimmed_1.fastq")
samples <- gsub("_trimmed_1.fastq.gz", "", samples)

# Conda Activate
conda_activate <- Sys.getenv("CONDA_ACTIVATE")

# Conda environment 
conda_env_wes <- paste0(Data_shared, Sys.getenv("CONDA_ENV_WES"))

# Mouse genome
genome <- paste0(Data_shared, Sys.getenv("GENOME"))

# Germline source
snps_indels_data <- paste0(Data_shared,Sys.getenv("SNPS_INDELS_DATA"))

# Interval list 
interval_list <- paste0(Data_shared, Sys.getenv("INTERVAL_LIST"))

# Bed file
Bed_file <-  paste0(Data_shared,Sys.getenv("BED_FILE"))

# Mapped sequences directory
Mapped_directory <-  paste0(Base_dir, Sys.getenv("MAPPED_DIRECTORY"))

# Duplicate metrics directory
Marked_dup_Directory <- paste0(Base_dir,Sys.getenv("MARKED_DUPLICATES_DIR"))

# Recall directory
Recall_directory <-  paste0(Base_dir,Sys.getenv("RECALL_DIR"))

# Jobs dir
jobs_dir <- paste0(Base_dir, Sys.getenv("JOBS_DIR"))

# Mapping QC directory
mapping_qc_dir <- paste0(Base_dir, Sys.getenv("MAPPED_QC_DIRECTORY"))

# Coverage file
Coverage_file <- paste0(Base_dir, Sys.getenv("COVERAGE_FILE"))

# Jobs requirements
cpu <- Sys.getenv("CPU_MAPPING")
memory <- Sys.getenv("MEMORY_MAPPING")
time <- Sys.getenv("TIME_MAPPING")
partition <- Sys.getenv("PARTITION_MAPPING")
mail <- Sys.getenv("MAIL")

################################################################################


################################################################################
#|  JOBS LOOP  
################################################################################

for(s in 1:length(samples)){

  # Name of the file (regarding the sample variable)
  filename <- paste(jobs_dir , "02_Mapping_", samples[s], ".sh", sep = "") 

  # Forward strand
  input_fastq_read1 <- paste(FATSQs_directory, samples[s], "_trimmed_1.fastq.gz", sep = "")

  # Reverse strand
  input_fastq_read2 <- paste(FATSQs_directory, samples[s], "_trimmed_2.fastq.gz", sep = "")
  
  # Output bam
  output_bam <- paste(Mapped_directory, samples[s], "_mapped.bam", sep = "") 

  # Output directory to store with added read groups
  output_rg_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.bam", sep = "")

  # Output directory to store sorted bam files
  output_rg_sorted_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.sorted.bam", sep = "") 

  # Output directory to store sorted bam files
  output_rg_sorted_dup_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.sorted.dup.bam", sep = "") 

  # Output directory to store sorted bam files
  output_rg_sorted_dup_bqsr_bam <- paste(Mapped_directory, samples[s], "_mapped.rg.sorted.dup.bqsr.bam", sep = "") 

  # Recal datafile
  recal_data_file <- paste(Recall_directory,"recal_data_",samples[s],".table", sep ="")

  cat(
    "#!/bin/bash",
    paste("#SBATCH --job-name=", samples[s], sep = ""),
    paste("#SBATCH -o logs/02_Mapping_", samples[s], ".out", sep = ""),
    paste("#SBATCH -e logs/02_Mapping_", samples[s], ".err", sep = ""),
    paste("#SBATCH --cpus-per-task=", cpu, sep = ""),
    paste("#SBATCH --mem=", memory, sep = ""),
    paste("#SBATCH --time=", time, sep = ""),
    paste("#SBATCH --partition=", partition, sep = ""),
    paste("#SBATCH --mail-user=",mail,sep = ""),
    "#SBATCH --mail-type=FAIL",

    "\n\n",
    paste("source ",conda_activate, sep = ""), 
    paste("conda activate ", conda_env_wes, sep = ""), 
    "\n\n",

    #"#---------------------------------- ALIGNMENT WITH BWA ----------------------------------",
    "\n",
    "# STEP 1 and 2: BWA TO ALIGN FASTQ FILES WITHOUT SAVING SAM FILES",
    paste("bwa mem -t ",cpu, " ", genome," ", input_fastq_read1, " ", input_fastq_read2, " | samtools view -@",cpu," -S -b -o ", output_bam, sep = ""),

    "\n",
    "# STEP 3: Add read groups",
    paste("picard AddOrReplaceReadGroups I=", output_bam," O=", output_rg_bam, " RGID=Sample_",samples[s], " RGLB=Library1 RGPL=illumina RGPU=Unit1 RGSM=Sample_",samples[s], sep = ""),

    "\n",
    "# STEP 4: Sort BAM USING SAMTOOLS",
    paste("samtools sort -@",cpu," ", output_rg_bam, " -o ", output_rg_sorted_bam, sep = ""),
    paste("samtools index", output_rg_sorted_bam, sep = " "),

    "\n\n",
   
    "#------------------------------------ MARK DUPLICATES -----------------------------------",
    "\n",

    paste("picard MarkDuplicates INPUT=", output_rg_sorted_bam," OUTPUT=", output_rg_sorted_dup_bam, " METRICS_FILE=", Marked_dup_Directory , "marked_dup_metrics_",samples[s], ".txt CREATE_INDEX=TRUE REMOVE_DUPLICATES=TRUE", sep = ""),

    "\n\n",

    "#--------------------------- BASE QUALITY SCORE RECALIBRATION ---------------------------",
    "\n",

    "# Step 1: Build the model",
    paste("gatk BaseRecalibrator -I", output_rg_sorted_dup_bam, 
    "-R", genome, 
    "-L", interval_list, 
    "--known-sites", snps_indels_data, 
    "-O", recal_data_file, sep = " "),
    "\n",

    "# Step 2: Apply the model to adjust the base quality scores",
    paste("gatk ApplyBQSR -I", output_rg_sorted_dup_bam, 
    "-R", genome, 
    "-L", interval_list, 
    "--bqsr-recal-file", recal_data_file,
    "-O", output_rg_sorted_dup_bqsr_bam, sep = " "),
    "\n",
    "# Step 3: Indexing the final results",
    paste("samtools index", output_rg_sorted_dup_bqsr_bam, sep = " "),

    "\n\n",

    "#--------------------------------- MAPPING QC METRICS ----------------------------------",
    "\n",
    
    "# Step1: Mapping QC using qualimap",
    paste("mkdir ", mapping_qc_dir,"qualimap_",samples[s], sep = ""),
    
    paste("samtools flagstat ", output_rg_sorted_dup_bqsr_bam, " > ", mapping_qc_dir,samples[s],"_samtools_qc.txt", sep =""),

    #paste("qualimap bamqc --paired --skip-duplicated --sorted -bam ", output_rg_sorted_dup_bqsr_bam, " -outdir ", mapping_qc_dir,"qualimap_",samples[s],"/", sep = ""),
    
    #"\n",
    #paste("mean_depth=$(bedtools coverage -a ",Bed_file," -b ", output_rg_sorted_dup_bqsr_bam, " -d | awk '{sum+=$5} END {if (NR>0) print sum/NR; else print '0'}')", sep = ""),
    #paste("echo -e", samples[s],"\t${mean_depth} >> ",samples[s], "_",Coverage_file, sep = ""),

    file = filename, sep = "\n", append = FALSE)

  system(paste("sbatch", filename, sep = " "))

}