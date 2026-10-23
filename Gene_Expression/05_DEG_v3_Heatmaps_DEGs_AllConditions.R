################################################################################
#                 DEGs HEATMAP/PCA IN ALL SAMPLES 
################################################################################

# Date: Thu Sep 25 12:35:33 2025
# Author: Maria Ponce



################################################################################
#                         PARAMETERS SELECTION 
################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "project"

# Pathway to the folders and files
path <- "path"

# Date of the log file 5_DEG_qc_XXXX.log
analysis_ID <- "analysis_ID" # Combatseq; ~ Condition; p = 0.05

# Select analysis/method performed
analysis <- "DESeq2"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
source(paste(path, project, "libraries_degs.R", sep = ""))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Set working directories ##

# General output directory
dir_out <- paste(path, project, sep = "")

# Output folder
dir_output <- paste(dir_out, "/05_DEG_ANALYSIS/", analysis_ID, sep = "")

# Files in 
dir_infiles <- paste(dir_output, "/Results/", sep = "")


# Load log file 
logfile <- read.table(paste(dir_out, "/log/5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)

# Experimental condition
trt <- logfile$condition

# Contrast levels
lvl_ord <- unlist(str_split(logfile$condition_order, pattern = ","))

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

### Pre-processing cutoffs

# Outliers 
# Remove the samples considered outliers
if(is.na(logfile$Outliers)){
  outliers <- NULL
} else {
  outliers <-  unlist(str_split(logfile$Outliers, ","))
  project <- paste(project, paste(outliers, collapse = "_"), sep = "_")}


# Variance transformation method used 
vsd_type <- logfile$Variance


### Threshold criteria 
fdr_cutoff <- logfile$fdr_cutoff
lfc_cutoff <- logfile$lfc_cutoff

# Threshold label
threshold <- paste("padj_", fdr_cutoff, "_log2FC_", round(lfc_cutoff, 2), sep ="")

# Color list
cols <- unlist(str_split(logfile$colortrt, pattern = ","))
names(cols) <- lvl_ord


################################################################################
#                         Load data  
################################################################################

# Load metadata
metadata <- read.delim(paste(dir_infiles, "Metadata_", project, "_", analysis_ID, ".txt", sep = ""), header = TRUE, sep = " ")

# Load DEGs results
data <- read.table(paste(dir_infiles, analysis, "_", project, ";Selected_", vsd_type, "blindFALSE_", threshold, "_", analysis_ID, ".txt", sep =""), dec = ".", header = TRUE)
print(dim(data))


###############################################################
#                           PLOTS
###############################################################

for (i in 1:length(contrast)) {
  
  # Contrast 
  name <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
  print(name)
  
  # Comparison folder
  dir_fig <- paste(dir_output, name, analysis, sep = "/")
  setwd(dir_fig)
  
  # Order levels according to the comparison
  ordered_lvl <- c(contrast[[i]][3], contrast[[i]][2])
  ordered_lvl <- c(ordered_lvl, lvl_ord[-which(lvl_ord %in% ordered_lvl)])
  
  
  # Step 1: Select comparison of interest
  # ---------------------------------------------------------------
  df <- data[which(data$Comparison == name), ]
  df_sel <- df
  
  # Step 2: Select comparison of interest
  # ---------------------------------------------------------------
  metadata <- metadata[which(metadata[[trt]] %in% ordered_lvl),]
  metadata[[trt]] <- factor(metadata[[trt]], levels = ordered_lvl, ordered = TRUE)
  metadata <- metadata[order(metadata[[trt]], metadata$Sample),]
  metadata$Color <- NA
  
  
  for(i in 1:length(ordered_lvl)){
    sel_level <- ordered_lvl[i]
    sel_color <- cols[which(names(cols) %in% sel_level)]
    metadata$Color[which(metadata[[trt]] %in% sel_level)] <- sel_color
  }
  
  print(head(metadata))
  
  # Sample information for the heatmap with colours
  sample_col <- data.frame(Condition = metadata[[trt]])
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
           clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/Heatmap_zscore_AllConditions_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), 
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
           clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/Heatmap_zscore_AllConditions_colors_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), 
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
           clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/Heatmap_zscore_AllConditions_woclustercols_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), 
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
           clustering_callback = callback,     # Sorting dendrogram
           fontsize = 6,
           filename = paste(dir_fig, "/Heatmap_zscore_AllConditions_woclustercols_colors_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), 
           height = 4, width = 4)
  
  
  
  # PCA plots
  #
  # To perform the PCA analysis we need a matrix with the PSI values per each 
  # sample. The samples must be place in the rows and the events in the columns. 
  # This matrix is centered and scale by using prcomp().
  #
  # A data frame with the original data and the metadata must be used to plot.
  # This data is used to increase interpretation power.
  #
  # The principal components which explains the variability of the data should
  # at least sum 65% of the variance of the data.
  
  m_t <- t(m)
  print(dim(m_t))
  m_pca <- prcomp(m_t, scale. = TRUE)
  
  var_perc <- get_eigenvalue(m_pca)[,2]
  
  # Scree plot
  # Percentage of variances explained by each principal component
  pca_scree <- fviz_eig(m_pca, choice = "variance", addlabels = TRUE, ggtheme = theme_classic(), main = "Screeplot")
  
  # PC1 vs PC2
  pca_1vs2 <- fviz_pca_ind(m_pca, axes = c(1, 2),
                           geom.ind = c("point", "text"), 
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]],  fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           title = "") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC2 (", round(var_perc[2],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # PC1 vs PC3
  pca_1vs3 <- fviz_pca_ind(m_pca, axes = c(1, 3),
                           geom.ind = c("point", "text"),  
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE, 
                           col.ind = metadata[[trt]], fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC3 (", round(var_perc[3],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # PC1 vs PC4
  pca_1vs4 <- fviz_pca_ind(m_pca, axes = c(1, 4),
                           geom.ind = c("point", "text"),
                           pointshape = 21, labelsize = 4, repel = TRUE, mean.point = FALSE,
                           col.ind = metadata[[trt]], fill.ind = metadata[[trt]],
                           addEllipses = FALSE, ellipse.level = 0.95,
                           legend.title = "Treatment", title = "") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = paste("PC1 (", round(var_perc[1],2),"% of variance)", sep = ""), y = paste("PC4 (", round(var_perc[4],2),"% of variance)", sep = ""))+
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text = element_text(size = 6),
          panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
  
  # Save pcas
  ggsave(filename = paste("PCA_params_AllConditions_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), plot = pca_scree, path = dir_fig, height = 4, width = 4, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs2)), "_AllConditions_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), plot = pca_1vs2, path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs3)), "_AllConditions_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), plot = pca_1vs3,  path = dir_fig, height = 5, width = 6, bg = "white")
  ggsave(filename = paste(deparse(substitute(pca_1vs4)), "_AllConditions_", name, "_n_", nrow(m), "_", project, "_", analysis_ID, ".pdf", sep = ""), plot = pca_1vs4, path = dir_fig, height = 5, width = 6, bg = "white")
  
  
}


