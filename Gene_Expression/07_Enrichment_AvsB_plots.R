################################################################################
#                 ENRICHMENT ORA DOTPLOTS
################################################################################

# Date: 04/09/2025
# Author: Maria Ponce

# Summary
# ---------------------
# 
# Generate top 10 dotplots based on the recall/precision per each sheet of the 
# enrichment analysis. Considered that each source is considered independent.



################################################################################
#                             LOAD LIBRARIES 
################################################################################


library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(openxlsx)


################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
project <- "project"

# Pathway to the folders and files
path <- "path"


# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID <- "analysis_ID" # Sample_info QC PDF; ~ Replicate + Condition

# Select analysis/method performed
analysis <- "DESeq2"

# Number of top pathways 
n <- 10
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


################################################################################
#                         Create working directories
################################################################################

# Create a folder for the enrichment analysis results. This folder contains 
# several folders classified in Results and Figures. 

# General output directory
dir_out <- paste(path, project, sep = "")

# Directory libraries and log
dir_log <- paste(dir_out, "/log/", sep = "")

# Create enrichment directory
dir_output <- paste(dir_out, "/07_ENRICHMENT", sep='')

# Load output directory
dir_outfolder <- paste(dir_output, "/", analysis, "_", analysis_ID, sep='')

# Save files with comparison results separate
dir_res <- paste(dir_outfolder, "/Results", sep='')

# Figures folder
dir_fig <- paste(dir_outfolder ,"/Figures", sep='')
foldername <- paste("ORA_Top_", n, sep = "")
dir.create(file.path(dir_fig, foldername))
dir_fig <- paste(dir_fig ,"/", foldername, sep='')



# Load log file 
logfile <- read.table(paste(dir_log, "7_Enrichment_DESeq2_", analysis_ID, ".log", sep = ""), header = TRUE)

# Significance level
fdr_cutoff <- logfile$fdr_cutoff

# Log2 fold change threshold
lfc_cutoff <- round(logfile$lfc_cutoff, 2)

# Threshold label
theshold <- paste("padj_", fdr_cutoff, "_log2FC_", lfc_cutoff, sep = "")



################################################################################
#                           DOTPLOT FUNCTION  
################################################################################


dotplot_function <- function(data, variable, db, sheet){
  
  # Level order
  level_ord <- data$term_name
  
  # Figures labels
  labels <- ifelse(variable == "precision", "Enrichment Factor", "Gene Ratio")
  
  # Change term name
  if(sum((duplicated(level_ord) == TRUE))>0){
    # Change term name
    data$term_name <- paste(data$term_name, data$term_id, sep = "\n")
    # Level order
    level_ord <- data$term_name
  }
  
  # Dotplot
  ggplot(data = data, aes(y = factor(term_name, levels = rev(level_ord)), x = !!as.name(variable), size = !!as.name(variable), color = p_value))+
    geom_point()+
    scale_color_gradient(name = "Adj p-value", low = "blue", high = "red")+
    scale_size(name =(labels))+
    labs(x = "", y = "", title = paste(sheet, ": ", db, " - ", variable, sep = ""))+ 
    theme_bw()+
    theme(axis.text.x = element_text(color = "black", size = 9), axis.text.y = element_text(color = "black", size = 9),
          legend.title = element_text(size = 9))
}



################################################################################
#                         Load data  
################################################################################


# File name 
filename <- paste(dir_res, list.files(dir_res), sep = "/")

# Sheet names
list_sheet <- excel_sheets(filename)

# Create a workbook
exc_rec <- createWorkbook()
exc_prec <- createWorkbook()

# Plots per sheet
for (sheet in list_sheet) {
  
  print(sheet)

  # Load data
  data <- read_xlsx(filename, sheet = sheet)
  
  # Database 
  database <- unique(data$source)
  
  # Create comparison folder
  dir.create(file.path(dir_fig, sheet))
  dir_sheet <- paste(dir_fig, sheet, sep = '/')
  
  # Data frame to save data
  recall <- c()
  precision <- c()
  
  
  # Loop for database
  for (db in database) {
    
    print(db)
    
    # Filter data
    df <- data[which(data$source == db),]
    print(dim(df))
    
    
    # Option 1: Sort by precision
    # ------------------------------------
    
    dt <- df[order(df$precision, df$p_value, decreasing = TRUE),]
    
    # Select top 15
    if(nrow(dt) >= n){dt <- dt[1:n,]}
    
    # Add data
    precision <- rbind(precision, dt)
    
    # Fig labels
    variable <- "precision"
    lbl <- paste(variable, sheet, db, analysis_ID, sep = "_")
    lbl <- gsub( ":", "", lbl) 
     
    # Dotplot 
    gpt <- dotplot_function(dt, variable, db, sheet)
    ggsave(filename = paste("00_", lbl, ".pdf", sep = ""), plot = gpt, path = dir_sheet, bg = "white", width = 8, height = 6)
    
    
    # Option 2: Sort by recall
    # ------------------------------------    
    
    dt <- df[order(df$recall, df$p_value, decreasing = TRUE),]
    
    # Select top 15
    if(nrow(dt) >= n){dt <- dt[1:n,]}
    
    # Add data
    recall <- rbind(recall, dt)
    
    # Fig labels
    variable <- "recall"
    lbl <- paste(variable, sheet, db, analysis_ID, sep = "_")
    lbl <- gsub( ":", "", lbl) 
    
    # Dotplot 
    gpt <- dotplot_function(dt, variable, db, sheet)
    ggsave(filename = paste("01_", lbl, ".pdf", sep = ""), plot = gpt, path = dir_sheet, bg = "white", width = 8, height = 6)
  }
  
  
  # Precision data
  addWorksheet(exc_prec, sheet)
  writeData(exc_prec, sheet, precision)
  
  # Recall data
  addWorksheet(exc_rec, sheet)
  writeData(exc_rec, sheet, recall)
  
}

# Save data
saveWorkbook(exc_prec, file = paste(dir_fig, "/Top_", n,"_ORA_precision_", project, "_", analysis_ID,".xlsx", sep = ""), overwrite = TRUE)
saveWorkbook(exc_rec, file = paste(dir_fig, "/Top_", n,"_ORA_recall_", project, "_", analysis_ID,".xlsx", sep = ""), overwrite = TRUE)
