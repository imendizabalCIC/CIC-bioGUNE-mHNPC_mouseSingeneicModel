################################################################################
#                                 GSEA
################################################################################

# Date: Wed Feb 19 15:37:11 2025
# Author: Maria Ponce

# Summary
# ---------------------
# 
# In this script, we prepare input data for GSEA. We have a pdf in which the format 
# is more detailed.
# 


################################################################################
#                             LOAD LIBRARIES 
################################################################################

library(fgsea)
library(openxlsx)
library(dplyr)
library(stringr)


################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "project"

# Pathway to the folders and files
path <- "path"


# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID <- "analysis_ID" # Sample_info QC PDF; ~ Replicate + Condition

# Select analysis/method performed
analysis <- "DESeq2"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# General output directory
dir_out <- paste(path, project, sep = "")


##############################################################################
#                         Create working directories
##############################################################################

# Create a folder for the enrichment analysis results. This folder contains 
# several folders classified in Results and Figures. 

# Create enrichment directory
dir.create(file.path(dir_out , "07_ENRICHMENT"), showWarnings = FALSE)
dir_output <- paste(dir_out, "/07_ENRICHMENT", sep='')

# Analysis ID data
dir.create(file.path(dir_output, paste(analysis, analysis_ID, sep = "_")), showWarnings = FALSE)
dir_outfolder <- paste(dir_output, "/", analysis, "_", analysis_ID, sep='')

# Create GSEA directory 
dir.create(file.path(dir_outfolder, "GSEA"), showWarnings = FALSE)
dir_outfolder <- paste(dir_outfolder, "/GSEA", sep='')
setwd(dir_outfolder)



################################################################################
#                               SPLICING 
################################################################################

if(analysis == "Splicing"){
  
  # Load log file 
  logfile <- read.table(paste(dir_wd, "/log/4_AS_qc_", logdate, ".log", sep = ""), header = TRUE)
  
  # Input directory 
  dir_trajectories <- paste(dir_out, "/06_TRAJECTORIES/", analysis, "/", clustering, "/Results", sep = "")
  
  # Output directory
  dir_infiles <- paste(dir_out, "/04_SPLICING_ANALYSIS/Results", sep = "")
  
  # Select the identification for the analysis
  label_gen <- "GENE"
  
  # Outliers 
  # Remove the samples considered outliers
  if(is.na(logfile$Outliers)){outliers <- NULL} else {outliers <-  paste(logfile$Outliers, collapse = ",")}
  
  # Experimental condition
  # Choose only one condition per script
  # The name should be the same as in the metadata file or sample information file
  trt <- logfile$condition
  
  # Contrast levels
  # The order is important because it will be used to order the data, as well as, 
  # to create the contrast matrix, reference of the order, plot data, ...
  # 
  # The first must be the reference level
  lvl_ord <- as.numeric(unlist(str_split(logfile$condition_order, pattern = ",")))
  
  # Contrast
  contrast <- unlist(str_split(logfile$contrast, ","))
  contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))
  
  
  # Threshold criteria 
  # Significance level
  # Default. pvalue < 0.05 
  pval_cutoff <- logfile$pval_cutoff
  # Log2 fold change threshold
  # Default. log2(1.5)
  DPSI_cutoff <- logfile$DPSI_cutoff
  
  # Threshold label
  threshold <- threshold <- paste("DPSI_", DPSI_cutoff ,"_pval_", pval_cutoff, sep = "")
  
  # File name
  if(clustering == "FCM"){
    load_file <- paste(clustering, "_", analysis, "_", project, ";", theshold, "_m_", m,"_c_", c, sep = "")
  } else {load_file <- paste(clustering, "_", analysis, "_", project, ";", theshold, "_c_", c, sep = "")}
  
  
  ################################################################################
  #                                 DEGS 
  ################################################################################
  
} else {
  # Load log file 
  logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)
  
  # Output directory
  dir_infiles <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID,"/Results", sep = "")
  
  # Select the identification for the analysis
  label_gen <- "Ensembl"
  
  # Condition
  trt <- logfile$condition
  
  # Significance level
  fdr_cutoff <- logfile$fdr_cutoff
  
  # Log2 fold change threshold
  lfc_cutoff <- round(logfile$lfc_cutoff, 2)
  
  # Organism
  organism <- logfile$Organism
  
  # Threshold label
  theshold <- paste("padj_", fdr_cutoff, "_log2FC_", lfc_cutoff, sep = "")
  
  # Variance stabilization
  md <- logfile$Variance
  
  # Contrast
  contrast <- unlist(str_split(logfile$contrast, ","))
  contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))
  
  ### Pre-processing cutoffs
  
  # Outliers 
  # Remove the samples considered outliers
  if(is.na(logfile$Outliers)){
    outliers <- NULL
  } else {
    outliers <-  paste(logfile$Outliers, collapse = ",")
    project <- paste(project, paste(outliers, collapse = "_"), sep = "_")}
  
  
  
  ################################################################################
  #                         Load data  
  ################################################################################
  
  # File name 
  load_file <- paste(analysis, "_", project, ";All_", md, "blindFALSE_", theshold, sep = "")
  
  # Load universe
  bg <- read.table(paste(dir_infiles, "/", load_file, "_", analysis_ID, ".txt", sep = ""), header = TRUE)
  
}

# Load metadata
sample_info <- read.table(paste(dir_infiles, "/", "Metadata_", project, "_", analysis_ID,".txt", sep = ""), header = TRUE)



################################################################################
#                               PROCESS DATA
################################################################################


# Assembly
if(organism == "Human"){assembly <- "hsapiens"}else{assembly <- "mmusculus"}


for (i in 1:length(contrast)) {
  
  # Contrast 
  name <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
  print(name)
  
  # Levels 
  lvl_ord <- c(contrast[[i]][3], contrast[[i]][2])
  
  
  if(length(which(bg$Comparison == name)) != 0){
  
    # 1. Prepare phenotype data 
    # --------------------------------------
    
    # Metadata 
    metadata <- sample_info[which(sample_info$Condition %in% lvl_ord), ]
    
    # Save data
    filename <- paste(name, "_phenotypedata_", analysis_ID, ".cls", sep = "");
    cat(paste(nrow(metadata), 2, 1, sep =" "),
        paste("#", paste(lvl_ord, collapse = " "), sep =" "),
        paste(metadata[[trt]], collapse = " "),
        file = filename, sep = "\n", append = FALSE)
    
    
    
    
    # 2. Prepare expression data 
    # --------------------------------------
    # Select the genes
    genes <- bg[which(bg$Comparison == name), c("Ensembl", paste("Norm", metadata$Sample, sep = "_"))]
    colnames(genes) <- gsub("Norm_", "", colnames(genes))
    
    if(identical(colnames(genes), c("Ensembl", metadata$Sample) == FALSE)){
      genes <- genes[, match(colnames(genes), c("Ensembl", metadata$Sample))]
    }
    
  
    # Process expression data
    colnames(genes)[1] <- "NAME"
    exprss_data <- genes %>% mutate(DESCRIPTION = "na") %>% select(NAME, DESCRIPTION, everything())
    
    # Save data
    write.table(exprss_data, paste(name, "_expressiondata_", analysis_ID, ".txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
    
  } else {sprintf("No DEG was found in comparison %s", name)}
  
}

