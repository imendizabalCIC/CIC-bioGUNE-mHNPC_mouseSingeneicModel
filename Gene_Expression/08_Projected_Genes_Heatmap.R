################################################################################
#                              HEATMAP PLOT
################################################################################

# Date: Wed Jul 16 11:30:03 2025
# Author: Maria Ponce

# Summary
# ---------------------

# This script creates the heatmap plots with all the significant genes.

# Please use the same size as indicated bellow
# height = 4 and width = 4
# 
# ### Access R ###
# R
# # Load R libraries 
# .libPaths("LibPaths")


################################################################
#                     LOAD LIBRARIES 
################################################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(require(openxlsx))
suppressPackageStartupMessages(library(pheatmap))  
suppressPackageStartupMessages(library(dendsort))  



################################################################
#                       SELECT PARAMETERS
################################################################

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PARAMETERS MUST BE SPECIFY PER EACH PROJECT
# Use to load/save the data of each specific project
# 
path <- "path"

# Project name 
project <- "project"

# Analysis ID
analysis_ID <- "analysis_ID" 

# DESeq2/DESeq2_nofilter 
analysis <- "DESeq2"

# Contrast
comp <- "VehiclevsDHT"
# Contrast level order
comp_lvl <- c("Vehicle", "DHT", "ENZA", "ENZA_DHT")
# Colors
cols <- c("#D1BDA2FF", "#74070E", "#B0A6DF", "#4E6813")
# # Contrast level order
# comp_lvl <- c("Vehicle", "DHT", "ENZA_DHT")
# # Colors
# cols <- c("#D1BDA2FF", "#74070E", "#4E6813")
names(cols) <- comp_lvl

# remove variance=0 rows
novar.fx <- function(x) { return(var(x,na.rm=T))}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sorting the dendrogram
# Modify parameter clustering_callback
# 
# The columns must followed this sorting:
# The control/vehicle samples must be localized on the left size of the heatmap
# The Condition/SAM486A samples must be localized on the right side of the heatmap

callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")} # Default in the pheatmap
# callback = function(hc, ...){dendsort(hc, isReverse = TRUE, type = "average")}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



################################################################################
#                     SET WORKING DIRECTORY
################################################################################


# Output directory
dir_out <- paste(path, project, "/", sep = "")

# Logfile directory
dir_log <- paste(dir_out, "log/", sep = "")

# Input directory
dir_input <- paste(dir_out, "05_DEG_ANALYSIS/", analysis_ID, "/Results/", sep = "")

# Splicing paper figures directory
dir.create(file.path(dir_out, "08_Projected_Genes"))
dir_folder <- paste(dir_out, "08_Projected_Genes", sep = "")
dir.create(file.path(dir_folder, analysis_ID))
dir_folder <- paste(dir_folder, analysis_ID, sep = "/")



################################################################
#                       LOGFILE
################################################################

# Load logfile
logfile <- read.table(paste(dir_log ,"5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)

# Variance stabilization method
vsd_type <- logfile$Variance

### Threshold criteria 

## Significance level
# Default. pvalue < 0.05 
fdr_cutoff <- logfile$fdr_cutoff

## Log2 fold change threshold
# Default. log2(1.5)
lfc_cutoff <- logfile$lfc_cutoff

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

## Threshold label
threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")




################################################################
#                       LOAD DATA
################################################################

# Load metadata
metadata <- read.delim(paste(dir_input, "Metadata_", project, "_", analysis_ID, ".txt", sep = ""), header = TRUE, sep = " ")

# Load Matt results
data <- read.table(paste(dir_input, analysis, "_", project, ";Selected_", vsd_type,"blindFALSE_", threshold, "_", analysis_ID, ".txt", sep =""), dec = ".", header = TRUE)
print(dim(data))


###############################################################
#                           PLOTS
###############################################################

for (i in 1:length(contrast)) {
  
  # Contrast 
  name <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
  print(name)
  
  # Comparison folder
  dir.create(file.path(dir_folder, name))
  dir_folder_comp <- paste(dir_folder, name, sep = "/")
  
  # PCA folder
  dir.create(file.path(dir_folder_comp, "Heatmap"))
  dir_fig <- paste(dir_folder_comp, "Heatmap", sep = "/")
  setwd(dir_fig)
  
  
  # Order levels according to the comparison
  ordered_lvl <- c(contrast[[i]][3], contrast[[i]][2])
  ordered_lvl <- c(ordered_lvl, comp_lvl[-which(comp_lvl %in% ordered_lvl)])
  
  
  # Step 1: Select comparison of interest
  # ---------------------------------------------------------------
  df <- data[which(data$Comparison == name), ]
  df_sel <- df
  
  # Step 2: Select comparison of interest
  # ---------------------------------------------------------------
  metadata <- metadata[which(metadata$Condition %in% ordered_lvl),]
  metadata$Condition <- factor(metadata$Condition, levels = ordered_lvl, ordered = TRUE)
  metadata <- metadata[order(metadata$Condition, metadata$Sample),]
  metadata$Color <- NA
  
  
  for(i in 1:length(ordered_lvl)){
    sel_level <- ordered_lvl[i]
    sel_color <- cols[which(names(cols) %in% sel_level)]
    metadata$Color[which(metadata$Condition %in% sel_level)] <- sel_color
  }
  
  print(head(metadata))
  
  # Sample information for the heatmap with colours
  sample_col <- data.frame(Condition = metadata$Condition)
  rownames(sample_col) <- metadata$Sample
  
  
  ########### MATRIX ###########
  
  # Matrix with the PSI values 
  # 
  # This matrix will be used in different steps such as the heatmap, PCA 
  # or the correlation matrix. The information must be found:
  # - Rows as events 
  # - Columns as Samples
  m <- df_sel[,which(colnames(df_sel) %in% paste(paste(vsd_type, "_", sep = ""), metadata$Sample, sep = ""))]
  rownames(m) <- df_sel$Name
  colnames(m) <- gsub(paste(vsd_type, "_", sep = ""), "", colnames(m))
  m <- m[, match(metadata$Sample, colnames(m))]
  
  ###########  HEATMAP  ########### 
  # 
  # Clustering of all the differentially spliced events with different layers 
  # of information such as the Condition and the event type. 
  # The heatmap function needs a 
  #   - Matrix data (m): Samples as columns and events as rows 
  #   - Annotation events (sample_row): Event type corresponding to each event
  
  pheatmap(m,
           scale = "row",                    # Scale and center data 
           main = "",
           color = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),              # Color matrix
           border_color = NA,
           show_rownames = FALSE,
           show_colnames = FALSE,
           cluster_cols = TRUE,
           cluster_rows = TRUE, 
           treeheight_row = 0,                 # Do not plot the row dendrogram
           clustering_method = "ward.D2",      # Clustering method 
           # clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/00_Heatmap_VST_Zscore_", name, "_n_", nrow(m), "_", project,".pdf", sep = ""), 
           height = 4, width = 4)
  
  pheatmap(m,
           scale = "row",                    # Scale and center data 
           annotation_col = sample_col, 
           annotation_names_col = FALSE,
           annotation_colors = list(Condition = cols),
           main = "",
           color = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),              # Color matrix
           border_color = NA,
           show_rownames = FALSE,
           show_colnames = TRUE,
           cluster_cols = TRUE,
           cluster_rows = TRUE, 
           treeheight_row = 0,                 # Do not plot the row dendrogram
           clustering_method = "ward.D2",      # Clustering method 
           # clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/01_Heatmap_VST_Zscore_", name, "_n_", nrow(m), "_", project,".pdf", sep = ""), 
           height = 4, width = 4)
  
  pheatmap(m,
           scale = "row",                    # Scale and center data 
           main = "",
           color = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),              # Color matrix
           border_color = NA,
           show_rownames = FALSE,
           show_colnames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = TRUE, 
           treeheight_row = 0,                 # Do not plot the row dendrogram
           clustering_method = "ward.D2",      # Clustering method 
           # clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/00_Heatmap_VST_Zscore_woclustercols_", name, "_n_", nrow(m), "_", project,".pdf", sep = ""), 
           height = 4, width = 4)
  
  pheatmap(m,
           scale = "row",                    # Scale and center data 
           annotation_col = sample_col, 
           annotation_names_col = FALSE,
           annotation_colors = list(Condition = cols),
           main = "",
           color = rev(colorRampPalette(c("red4", "snow1", "royalblue4"))(50)),              # Color matrix
           border_color = NA,
           show_rownames = FALSE,
           show_colnames = TRUE,
           cluster_cols = FALSE,
           cluster_rows = TRUE, 
           treeheight_row = 0,                 # Do not plot the row dendrogram
           clustering_method = "ward.D2",      # Clustering method 
           # clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/01_Heatmap_VST_Zscore_woclustercols_", name, "_n_", nrow(m), "_", project,".pdf", sep = ""), 
           height = 4, width = 4)
  
}















