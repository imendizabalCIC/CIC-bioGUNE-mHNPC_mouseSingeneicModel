################################################################################
#                          PLOT HALLMARK GENES
################################################################################

# Date: Tue Aug 12 11:59:41 2025
# Author: Maria Ponce

# Summary
# ---------------------
# 
# From GSEA data, select a hallmark of interest to plot all the genes in a heatmap 
# together across the comparisons and plot them separately.

################################################################
#                       SELECT PARAMETERS
################################################################

# Project name 
project <- "project"

# Path 
path <- "path"

# DESeq2/DESeq2_nofilter 
analysis <- "DESeq2"

# Select hallmark
hallmark <- "HALLMARK_ANDROGEN_RESPONSE"

# Analysis ID
analysis_ID <- "20250924145419" #
# GSEA folder of interest
gsea_folder <- "VehvsDHT.Gsea.1758725888984"
# Comparison 
comparison <- "VehiclevsDHT"
lvl_ord <- c("Vehicle", "DHT")


# Load libraries
source(paste(path, project, "Libraries_degs.R", sep = ""))


################################################################################
#                   SET WORKING DIRECTORIES
################################################################################

# Output directory
dir_out <- paste(path, project, "/", sep = "")

# Log directory
dir_log <- paste(path, project, "/log/", sep = "")

# Input directories  
dir_infiles <- paste(dir_out, "05_DEG_ANALYSIS/", analysis_ID, "/Results/", sep = "")

# Figure directory
dir_output <- paste(dir_out, "07_ENRICHMENT/", analysis, "_", analysis_ID, "/Figures/", sep = "")

# Figure folder 
dir.create(file.path(dir_output, hallmark))
dir_fig <- paste(dir_output, hallmark, sep = "")
setwd(dir_fig)


################################################################################
#                             DATA 
################################################################################

# Filename
filename <- list.files(dir_infiles, pattern = paste("DESeq2_", project, ";All_VSTblindFALSE", sep = ""))[1]

# Load data
data <- read.table(paste(dir_infiles, filename, sep = ""), header = TRUE)
print(dim(data))
print(head(data))

# Load metadata
metadata <- read.table(paste(dir_infiles, "Metadata_", project, "_", analysis_ID, ".txt", sep = ""), header = TRUE)

# Load log data
logfile <- read.table(paste(dir_log ,"5_DEG_qc_", analysis_ID, ".log", sep = ""), header = TRUE)


################################################################################
#                             HALLMARKS 
################################################################################

# Load hallmarks .gmt

# File format
# 
# Hallmark1 Description Gen1 Gen2 Gen3 
# Hallmark2 Description Gen1 Gen2 Gen3
# Hallmark3 Description Gen1 Gen2 Gen3

# Load file
gmtfile <- readLines(paste(dir_out, "07_ENRICHMENT/", analysis, "_", analysis_ID, "/GSEA/", gsea_folder, "/edb/gene_sets.gmt", sep = ""))
# [1] "HALLMARK_CHOLESTEROL_HOMEOSTASIS\thttps://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_CHOLESTEROL_HOMEOSTASIS\tAbca2\tAcat2\tAcss2\tActg1\tAlcam\tAldoc\tAntxr2\tAnxa5\tAtf3\tAtf5\tAtxn2\tCd9\tChka\tClu\tCpeb2\tCtnnb1\tCxcl16\tCyp51\tDhcr7\tEbp\tEch1\tErrfi1\tEthe1\tFabp5\tFads2\tFasn\tFbxo6\tFdft1\tFdps\tGldc\tGnai1\tGpx8\tGstm7\tGusb\tHmgcr\tHmgcs1\tHsd17b7\tIdi1\tJag1\tLdlr\tLgals3\tLgmn\tLss\tMal2\tMvd\tMvk\tNfil3\tNiban1\tNsdhl\tPcyt2\tPdk3\tPlaur\tPlscr1\tPmvk\tPnrc1\tPparg\tS100a11\tSc5d\tSema3b\tSqle\tSrebf2\tStard4\tStx5a\tTm7sf2\tTmem97\tTnfrsf12a\tTrib3\tTrp53inp1"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
# [2] "HALLMARK_MTORC1_SIGNALING\thttps://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MTORC1_SIGNALING\tAbcf2\tAcaca\tAcly\tAcsl3\tActr2\tActr3\tAdd3\tAdipor2\tAk4\tAldoa\tArpc5l\tAsns\tAtp2a2\tAtp5mc1\tAtp6v1d\tAurka\tBcat1\tBhlhe40\tBtg2\tBub1\tCacybp\tCalr\tCanx\tCcnf\tCcng1\tCct6a\tCd9\tCdc25a\tCdkn1a\tCfp\tCops5\tCoro1a\tCth\tCtsc\tCyb5b\tCyp51\tDapp1\tDdit3\tDdit4\tDdx39a\tDhcr24\tDhcr7\tDhfr\tEbp\tEdem1\tEef1e1\tEgln3\tEif2s2\tElovl5\tElovl6\tEno1b\tEprs1\tEro1a\tEtf1\tFads1\tFads2\tFdxr\tFkbp2\tG6pdx\tGapdh-ps15\tGbe1\tGclc\tGga2\tGla\tGlrx\tGmps\tGot1\tGpi1\tGsk3b\tGsr\tGtf2h1\tHk2\tHmbs\tHmgcr\tHmgcs1\tHprt1\tHsp90b1\tHspa4\tHspa5\tHspa9\tHspd1\tHspe1\tIdh1\tIdi1\tIfi30\tIfrd1\tIgfbp5\tImmt\tInsig1\tItgb2\tLdha\tLdlr\tLgmn\tLta4h\tM6pr\tMap2k3\tMcm2\tMcm4\tMe1\tMllt11\tMthfd2\tMthfd2l\tNampt\tNfil3\tNfkbib\tNfyc\tNherf1\tNiban1\tNmt1\tNufip1\tNup205\tNupr1\tP4ha1\tPdap1\tPdk1\tPfkl\tPgk1\tPgm1\tPhgdh\tPik3r3\tPitpnb\tPlk1\tPlod2\tPno1\tPnp\tPolr3g\tPpa1\tPpia\tPpp1r15a\tPrdx1\tPsat1\tPsma3\tPsma4\tPsmb5\tPsmc2\tPsmc4\tPsmc6\tPsmd12\tPsmd13\tPsmd14\tPsme3\tPsmg1\tPsph\tQdpr\tRab1a\tRdh11\tRit1\tRpa1\tRpn1\tRrm2\tRrp9\tSc5d\tSdf2l1\tSec11a\tSerp1\tSerpinh1\tShmt2\tSkap2\tSla\tSlc1a4\tSlc1a5\tSlc2a1\tSlc2a3\tSlc37a4\tSlc6a6\tSlc7a11\tSlc7a5\tSord\tSqle\tSqstm1\tSrd5a1\tSsr1\tStard4\tStc1\tStip1\tSytl2\tTbk1\tTcea1\tTes\tTfrc\tTm7sf2\tTmem97\tTomm40\tTpi1\tTrib3\tTuba4a\tTubg1\tTxnrd1\tUbe2d3\tUchl5\tUfm1\tUng\tUso1\tVldlr\tWars1\tXbp1\tYkt6"


## Process gmt file

# Generate a list per each line a separeta each element based on the separation mark ("\t")
gmt_list <- strsplit(gmtfile, "\t")

# Convert to a named list
# - name = gene set name
# - value = vector of genes
# 
# Do not consider the first two elements, which correspond to the hallmark and the descriptions
hallmark_list <- lapply(gmt_list, function(x) x[3:length(x)])
names(hallmark_list) <- sapply(gmt_list, function(x) x[1])

## Select hallmark of interest
genelist <- hallmark_list[[hallmark]]

# save gene list
writeLines(genelist, paste(hallmark, ".txt", sep =""))


################################################################################
#                             SET PARAMETER 
################################################################################


# Experimental condition
trt <- logfile$condition

# Contrast levels
lvl_ord <- unlist(str_split(logfile$condition_order, pattern = ","))

# Comparison vector
contrasts <- unique(data$Comparison)

# Threshold label
threshold <- paste("padj_", logfile$fdr_cutoff, "_log2FC_", round(logfile$lfc_cutoff, 2), sep ="")

# Color list
color_list <- list(trt = unlist(str_split(logfile$colortrt, pattern = ",")), 
                   Heatmap = unlist(str_split(logfile$colorheat, pattern = ",")),
                   Shared = unlist(str_split(logfile$colorsh, pattern = ",")))
names(color_list) <- c(trt, "Heatmap", "Shared")
names(color_list[[trt]]) <- lvl_ord

# Plot label 
label <- paste(hallmark, project, analysis_ID, sep = "_")


################################################################################
#                           PROCESS DATA 
################################################################################

## Process data 

# Remove duplicated genes across the comparison
df <- data[-which(duplicated(data$Name)),]
print(dim(df))

# Select genes from the hallmark
hm_df <- df[which(df$Symbol %in% genelist),]
dim(hm_df)


## Process Heatmap data

# Matrix 
m <- hm_df[,which(colnames(df) %in% paste("VST_", metadata$Sample, sep = ""))]
colnames(m) <- gsub("VST_", "", colnames(m))
rownames(m) <- hm_df$Symbol

# Sample information for the heatmap with colors
sample_col <- data.frame(Condition = metadata[[trt]])
rownames(sample_col) <- metadata$Sample


sample_col <- data.frame(Condition = metadata$Condition[which(metadata$Condition %in% lvl_ord)])
rownames(sample_col) <- metadata$Sample[which(metadata$Condition %in% lvl_ord)]
m <- m[, rownames(sample_col)]

################################################################################
#                            PLOT
################################################################################


# Sorting the dendrogram
# Modify parameter clustering_callback
# 
# The columns must followed this sorting:
# The control/vehicle samples must be localized on the left size of the heatmap
# The Condition/SAM486A samples must be localized on the right side of the heatmap

callback = function(hc, ...){dendsort(hc, isReverse = FALSE, type = "average")} # Default in the pheatmap
# callback = function(hc, ...){dendsort(hc, isReverse = TRUE, type = "average")}


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
         color = color_list$Heatmap,              # Color matrix
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE, 
         treeheight_row = 0,                 # Do not plot the row dendrogram
         clustering_method = "ward.D2",      # Clustering method 
         clustering_callback = callback,     # Sorting dendrogram
         fontsize = 6,
         filename = paste(dir_fig, "/00_Heatmap_VST_Zscore_", comparison, "_", label, ".pdf", sep = ""),
         height = 4, width = 4)

pheatmap(m,
         scale = "row",                    # Scale and center data 
         annotation_col = sample_col, 
         annotation_names_col = FALSE,
         annotation_colors = color_list,
         main = "",
         color = color_list$Heatmap,              # Color matrix
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = TRUE, 
         treeheight_row = 0,                 # Do not plot the row dendrogram
         clustering_method = "ward.D2",      # Clustering method 
         clustering_callback = callback,     # Sorting dendrogram
         fontsize = 6,
         filename = paste(dir_fig, "/01_Heatmap_VST_Zscore_", comparison, "_", label, ".pdf", sep = ""),
         height = 4, width = 4)

pheatmap(m,
         scale = "row",                    # Scale and center data 
         main = "",
         color = color_list$Heatmap,              # Color matrix
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = TRUE, 
         treeheight_row = 0,                 # Do not plot the row dendrogram
         clustering_method = "ward.D2",      # Clustering method 
         clustering_callback = callback,     # Sorting dendrogram
         fontsize = 6,
         filename = paste(dir_fig, "/00_Heatmap_VST_Zscore_woclustercols_", comparison, "_", label, ".pdf", sep = ""), 
         height = 4, width = 4)

pheatmap(m,
         scale = "row",                    # Scale and center data 
         annotation_col = sample_col, 
         annotation_names_col = FALSE,
         annotation_colors = color_list,
         main = "",
         color = color_list$Heatmap,              # Color matrix
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = TRUE, 
         treeheight_row = 0,                 # Do not plot the row dendrogram
         clustering_method = "ward.D2",      # Clustering method 
         clustering_callback = callback,     # Sorting dendrogram
         fontsize = 6,
         filename = paste(dir_fig, "/01_Heatmap_VST_Zscore_woclustercols_", comparison, "_", label, ".pdf", sep = ""),
         height = 4, width = 4)


pheatmap(m,
         scale = "row",                    # Scale and center data 
         annotation_col = sample_col, 
         annotation_names_col = FALSE,
         annotation_colors = color_list,
         main = "",
         color = color_list$Heatmap,              # Color matrix
         border_color = NA,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = TRUE, 
         treeheight_row = 0,                 # Do not plot the row dendrogram
         clustering_method = "ward.D2",      # Clustering method 
         clustering_callback = callback,     # Sorting dendrogram
         fontsize_row = 2,         
         fontsize = 6,
         filename = paste(dir_fig, "/02_Heatmap_VST_Zscore_GeneNames_", comparison, "_", label, ".pdf", sep = ""),
         height = 4, width = 4)

pheatmap(m,
         scale = "row",                    # Scale and center data 
         annotation_col = sample_col, 
         annotation_names_col = FALSE,
         annotation_colors = color_list,
         main = "",
         color = color_list$Heatmap,              # Color matrix
         border_color = NA,
         show_rownames = TRUE, 
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = TRUE, 
         treeheight_row = 0,                 # Do not plot the row dendrogram
         clustering_method = "ward.D2",      # Clustering method 
         clustering_callback = callback,     # Sorting dendrogram
         fontsize_row = 3,         
         fontsize = 6,
         filename = paste(dir_fig, "/02_Heatmap_VST_Zscore_woclustercols_GeneNames_", comparison, "_", label, ".pdf", sep = ""),
         height = 4, width = 4)



################################################################################
#                            BOXPLOT
################################################################################

dir.create(file.path(dir_fig, "Violin"))
dir_fig2 <- paste(dir_fig, "/Violin/", sep = "")

for (i in 1:nrow(hm_df)) {
  
  # Select gene 
  gen_hm <- hm_df[i, ]
  gene <- gen_hm$Symbol
  
  print(gene)

  # Matrix 
  m <- gen_hm[, which(colnames(gen_hm) %in% paste("VST_", metadata$Sample, sep = ""))]
  colnames(m) <- gsub("VST_", "", colnames(m))

  # Stack matrix
  sm <- stack(m)
  sm$Condition <- rep(metadata[[trt]])
  sm$Color <- rep(color_list[[trt]], each = 4)
  
  # 
  ggplot(sm, aes(x = factor(Condition, levels = lvl_ord), y = values, fill = Color))+
    geom_violin(width = 0.9, alpha = 0.7, linewidth = 0.05)+
    geom_boxplot(width = 0.07, alpha = 0.7, linewidth = 0.05)+
    geom_point(size = 0.5)+
    scale_fill_identity()+
    labs(x = "", y = "Counts\n(VST)", title = gene)
  ggsave(filename = paste(dir_fig2, gene, "_Violin_", label,".pdf", sep = ""), height = 4, width = 4, bg = "white")
}







