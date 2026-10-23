################################################################################
#                        TOP GENE CANDIDATES
################################################################################

# Date: Tue Aug 12 10:23:12 2025
# Author: Maria Ponce

# Summary
# ---------------------
# 
# Generate a top candidate list for each comparison. I decided to choose the top 
# 50 candidates that were upregulated and downregulated in each comparison



################################################################################
#                              Set parameters
################################################################################


path <- "path"

# Project name 
project <- "project"

# Analysis ID
analysis_ID <- "analysis_ID" 

# DESeq2/DESeq2_nofilter 
analysis <- "DESeq2"

# Number of candidates
n <- 50


################################################################################
#                       SET WORKING DIRECTORIES
################################################################################

# Project directory
dir_out <- paste(path, project, "/", sep = "")

# Utils directory
dir_utils <- paste(dir_out, "utils/", sep = "")

# Input directories  
dir_infiles <- paste(dir_out, "05_DEG_ANALYSIS/", analysis_ID, "/Results/", sep = "")

# Output folder
dir.create(file.path(dir_out, "09_Top_Genes_Candidates"))
dir_output <- paste(dir_out, "09_Top_Genes_Candidates/", sep = "")
setwd(dir_output)


################################################################################
#                          LOAD LIBRARIES
################################################################################

source(paste(dir_utils, "libraries_degs.R", sep = ""))


################################################################################
#                             LOAD DATA
################################################################################

# Filename
filename <- list.files(dir_infiles, pattern = paste("DESeq2_", project, ";Selected_VSTblindFALSE", sep = ""))

# Load data
data <- read.table(paste(dir_infiles, filename, sep = ""), header = TRUE)
print(dim(data))
print(head(data))
print(table(data$Comparison))


################################################################################
#                            TOP CANDIDATE LIST  
################################################################################

# Contrasts
contrasts <- unique(data$Comparison)

# Results workbook
exc <- createWorkbook()


# Loop per comparison
for(comparison in contrasts){
  
  # Comparison
  print(comparison)
  
  # Select comparison DEGs
  df <- data[which(data$Comparison == comparison), ]
  print(dim(df))
  
  # Select upregulated/downregulated genes with a Symbol name
  up <- df[which(df$Direction == "Upregulated" & df$Symbol != ""),]
  down <- df[which(df$Direction == "Downregulated" & df$Symbol != ""),]
  
  print(nrow(up))
  print(nrow(down))
  
  # Sort by log2FC and adjusted p-value
  idx_u <- order(up$log2FC, up$padj, decreasing = TRUE)
  idx_d <- order(down$log2FC, down$padj, decreasing = FALSE)
  up <- up[idx_u,]
  down <- down[idx_d,]
  
  # Check upregulated 
  print(head(up))
  print(tail(up))
  
  # Check downregulated
  print(head(down))
  print(tail(down))
  
  # Select top 50 candidates upregulated/downregulated and remove comparison col
  sel_up <- up[1:n, -(1:2)]
  sel_down <- down[1:n, -(1:2)]

  
  # Save data
  addWorksheet(exc, sheetName = paste("Up", comparison, sep = "_"))
  addWorksheet(exc, sheetName = paste("Down", comparison, sep = "_"))
  writeData(exc, sel_up, sheet = paste("Up", comparison, sep = "_"))
  writeData(exc, sel_down, sheet = paste("Down", comparison, sep = "_"))
  
}

saveWorkbook(exc, paste("Top_", n, "_Candidates;", analysis, "_", project, "_", analysis_ID, ".xlsx", sep = ""), overwrite = TRUE)



