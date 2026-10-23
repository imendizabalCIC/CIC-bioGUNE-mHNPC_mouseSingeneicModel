##############################################################
#                    PCA PLOT SCRIPT
################################################################

# Date: 12/2024 
# Author: María Ponce

# This script creates the PCA plots with all the significant genes.

# Please use the same size as indicated bellow
# height = 4 and width = 4


################################################################
#                     LOAD LIBRARIES 
################################################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(require(openxlsx))
suppressPackageStartupMessages(library(scales))



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
names(cols) <- comp_lvl

# remove variance=0 rows
novar.fx <- function(x) { return(var(x,na.rm=T))}


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

# Comparison folder
dir.create(file.path(dir_folder, comp))
dir_folder <- paste(dir_folder, comp, sep = "/")

# PCA folder
dir.create(file.path(dir_folder, "PCA"))
dir_fig <- paste(dir_folder, "PCA", sep = "/")
setwd(dir_fig)


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


# Plot labels
name <- comp

# Step 1: Select comparison of interest
# ---------------------------------------------------------------
df <- data[which(data$Comparison == comp), ]
df_sel <- df



# Step 2: Select comparison of interest
# ---------------------------------------------------------------
metadata$Condition <- factor(metadata$Condition, levels = comp_lvl)
metadata$Color <- NA

for(i in 1:length(comp_lvl)){
  sel_level <- comp_lvl[i]
  sel_color <- cols[which(names(cols) %in% sel_level)]
  metadata$Color[which(metadata$Condition %in% sel_level)] <- sel_color
}

print(head(metadata))




########### MATRIX ###########

# Matrix with the PSI values 
# 
# This matrix will be used in different steps such as the heatmap, PCA 
# or the correlation matrix. The information must be found:
# - Rows as events 
# - Columns as Samples
m <- df_sel[,which(colnames(df_sel) %in% paste(vsd_type, "_", metadata$Sample, sep = ""))]
rownames(m) <- df_sel$NAME
colnames(m) <- gsub(paste(vsd_type, "_", sep = ""), "", colnames(m))
m <- m[, match(metadata$Sample, colnames(m))]
colnames(m) <- paste(metadata$Condition, metadata$Sample, sep = "_")


range_row <- function(row) {
  # Range gives the minimum and the maximum values found
  range_row <- range(row, na.rm = TRUE)
  # range is Max - Min
  diff(range_row)
}

row_ranges <- apply(m, 1, range_row)

if(length(which(row_ranges == 0)) >=1){
  m <- m[which(row_ranges == 0),]
} else {print("No row presented a range 0.")}


# novar <- apply(m,1,novar.fx)
# novar.rows <- c(which(novar == 0),which(is.na(novar)));
# 
#do PCA
# pc_data <- m[-c(novar.rows),]
pc <- prcomp(t(na.omit(m)), center = TRUE,scale = TRUE)
summary(pc)
pc.var1 <- 100*round(pc$sdev[1]^2/sum(pc$sdev^2),3);
pc.var2 <- 100*round(pc$sdev[2]^2/sum(pc$sdev^2),3); 

#make figure
col.v <- metadata$Color
text.x <- (paste("PC1 (",pc.var1,"% of variance)",sep=''));
text.y <- (paste("PC2 (",pc.var2,"% of variance)",sep=''));
#just plot to get the dimensions
plot(pc$x[,1],pc$x[,2],cex.lab=1,bty="l",xlab=text.x,ylab=text.y,pch=19,cex=2.5,col=alpha(col.v,0.5),main="",type="n")
coord <- par("usr")
#plot with extended space

pdf(paste("00_PC1_PC2_PCA_", name, "_n_", nrow(m), "_", project, ".pdf",sep=''),height=4,width=4)
plot(pc$x[,1],pc$x[,2],cex.axis=0.9,bty="l",xlab=text.x,ylab=text.y,pch=19,cex=2.2,col=alpha(col.v,0.5),main="",
     xlim=c(coord[1]+coord[1]*0.1,coord[2]+coord[2]*0.1),
     ylim=c(coord[3]+coord[3]*0.1,coord[4]+coord[4]*0.1))
dev.off();
#pdf with only legend
pdf(paste("00_PC1_PC2_PCA_", name, "_n_", nrow(m), "_", project, "_legend.pdf",sep=''),height=4,width=4)
plot(pc$x[,1],pc$x[,2],cex.axis=0.9,bty="l",xlab=text.x,ylab=text.y,pch=19,cex=2.2,col=alpha(col.v,0.5),main="",type="n")
legend(x = "center", 
       pch = 19,cex=1, col=alpha(cols,0.5),pt.cex=2.5,
       legend=comp_lvl,bty="n")  
dev.off()



