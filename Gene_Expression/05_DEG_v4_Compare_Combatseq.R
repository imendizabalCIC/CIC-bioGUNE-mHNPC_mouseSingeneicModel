################################################################################
#                 COMPARE DEGS ANALYSIS VERSION
################################################################################

# Date: Thu Sep 25 14:03:28 2025
# Author: Maria Ponce

# Summary
# ---------------------
# 
# This script aims to automatized the differentially expressed genes (DEGs)
# analysis based on the RNA-seq for the AC laboratory. In this script, we 
# 
# This script is a general analysis of the differences between analysis in the 
# same project in which we check DEGs among the comparisons by using Venn diagrams, 
# correlation plots and giving a table with the genes shared.


################################################################################
#                     LOAD LIBRARIES 
################################################################################


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))



################################################################
#                       SELECT PARAMETERS
################################################################

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PARAMETERS MUST BE SPECIFY PER EACH PROJECT
# Use to load/save the data of each specific project
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Project name
project <- "project"

# Pathway to the folders and files
path <- "path"

# Analysis ID
analysis_ID1 <- "analysis_ID1" # Sample_info QC PDF; ~ Replicate + Condition
analysis_ID2 <- "analysis_ID2" # Combatseq; ~ Condition; p = 0.05

# Select analysis/method performed
analysis <- "DESeq2"
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Figures labels
label1 <- "Raw Counts"
label2 <- "Combatseq Corrected Counts"


# DESeq2/DESeq2_nofilter 
analysis <- "DESeq2"

# Contrast level order
comp_lvl <- c("Vehicle", "DHT", "ENZA", "ENZA_DHT")

# Colors
cols <- c("#D1BDA2FF", "#74070E", "#B0A6DF", "#4E6813")
names(cols) <- comp_lvl

# Correation colors
cor_cols <- c("#74070E",  "#6699FF", "#4E6813")
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



################################################################################
#                       LOAD FUNCTIONS
################################################################################

# Check congruence 
check_congruence <- function(x){
  a <- unlist(strsplit(x, split = ", "))
  if (all(a ==a[1])){
    return("Yes")
  } else {return("No")}
}


################################################################################
#                       SET WORKING DIRECTORIES
################################################################################

# Project directory
dir_out <- paste(path, project, sep = "")

# Directory degs
dir_output <- paste(dir_out, "/05_DEG_ANALYSIS/", sep ="")

# Directory files in 
dir_analysis_ID1 <- paste(dir_output, analysis_ID1, "/Results/", sep ="")
dir_analysis_ID2 <- paste(dir_output, analysis_ID2, "/Results/", sep ="")

# Directory folder
dir.create(file.path(paste(dir_output, analysis_ID2, sep = ""), paste("Comparison", gsub(" ", "_", label1), analysis_ID1, gsub(" ", "_", label2), analysis_ID2, sep = "_")), showWarnings = F)
dir_output <- paste(dir_output, analysis_ID2, "/Comparison_", gsub(" ", "_", label1), "_", analysis_ID1, "_", gsub(" ", "_", label2), "_", analysis_ID2, sep = "")

# Directory results 
dir.create(file.path(dir_output, "Results"), showWarnings = F)
dir_res <- paste(dir_output, "/Results", sep = "")


################################################################################
#                       LOAD DATA
################################################################################


# Load metadata
metadata <- read.delim(paste(dir_analysis_ID1, "Metadata_", project, "_", analysis_ID1, ".txt", sep = ""), header = TRUE, sep = " ")

# Load data 
files <- list.files(dir_analysis_ID1, pattern = paste(analysis, "_", project, ";All_VSTblindFALSE", sep = ""))
filename1 <- files[grep(".txt", files)]
data1 <- read.table(paste(dir_analysis_ID1, filename1, sep =""), dec = ".", header = TRUE)
print(dim(data1))

files <- list.files(dir_analysis_ID2, pattern = paste(analysis, "_", project, ";All_VSTblindFALSE", sep = ""))
filename2 <- files[grep(".txt", files)]
data2 <- read.table(paste(dir_analysis_ID2, filename2, sep =""), dec = ".", header = TRUE)
print(dim(data2))


################################################################################
#                       Log file
################################################################################

# Load logfile
logfile <- read.table(paste(dir_out ,"/log/5_DEG_v2_", analysis, "_", analysis_ID1, ".log", sep = ""), header = TRUE)

# Contrast
contrast <- unlist(str_split(logfile$contrast, ","))
contrast <- split(contrast, rep(1:(length(contrast)/3), each = 3))

# Thresholds
fdr_cutoff <- logfile$fdr_cutoff
lfc_cutoff <- logfile$lfc_cutoff

# Figures and files general label
labs <- paste(project, gsub(" ", "_", label1), analysis_ID1, gsub(" ", "_", label2), analysis_ID2, sep = "_")

# ggplot2 theme 
theme_DEGs <- theme_bw()+ theme(axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), panel.grid = element_blank())
theme_set(theme_DEGs)



################################################################################
#                             PLOTS
################################################################################


# Transcriptome genes 
transcriptome <- list("Raw Counts" = unique(data1$Name),
                      "Combatseq Corrected Counts" = unique(data2$Name))

ggvenn(transcriptome, 
       digits = 2,                                    # Two decimals
       fill_color = as.vector(cols[c(1,4)]),       # Fill color
       fill_alpha = 0.7,                            # Transparency, default
       stroke_size = 0.25,                          # Line thickness
       set_name_size = 2.5,                           # Default 
       text_size = 2                              # Default
)
ggsave(filename = paste("/Venn_diagram_Transcriptome_", labs, ".pdf", sep = ""), path = dir_output, height = 4, width = 4, plot = last_plot())



# Loop for all the comparison

exc <- createWorkbook()

for(i in 1:length(contrast)){
  
  
  # Select comparison of interest
  comparison <- paste(contrast[[i]][3], "vs", contrast[[i]][2], sep = "")
  print(comparison)
  
  # Comparison levels
  lvl_ord <- c(contrast[[i]][2], contrast[[i]][3])
  
  # Color 
  color_list <- list(Condition = cols[which(names(cols) %in% lvl_ord)], 
                     Correlation = cor_cols)
  
  
  
  ##########################################################################
  #                       SET WORKING DIRECTORY
  ##########################################################################
  
  # Folder directory
  dir.create(file.path(dir_output, comparison), showWarnings = F)
  dir_fig <- paste(dir_output, comparison, sep = "/")
  
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # If the figures names are too long hte file will not be save
  # We change Vechicle to V to reduce the number of characters
  
  # Figure plot
  label_plot <- gsub("Vehicle", "V", paste(comparison, labs, sep = "_"))
  #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  ##########################################################################
  #                          PROCESS DATA 
  ##########################################################################
  
  # Select data
  df_filt1 <- data1[which(data1$Comparison == comparison),]
  df_filt2 <- data2[which(data2$Comparison == comparison),]
  
  df_filt1$Project <- label1
  df_filt2$Project <- label2
  
  filtered <- rbind(df_filt1, df_filt2)
  
  degs <- filtered[which(filtered$DEG == "YES"),]
  up <- degs[which(degs$Direction == "Upregulated"), ]
  down <- degs[which(degs$Direction == "Downregulated"), ]
  
  # Gene lists
  gene_list <- split(degs$Name, degs$Project)
  up_genes <- split(up$Name, up$Project)
  down_genes <- split(down$Name, down$Project)
  
  
  
  # DUPLICATED TABLES
  # ------------------------------------------------------------------------
  
  # Duplicated events
  dup_tab <- as.data.frame(table(degs$Name))
  dup_event <- dup_tab[dup_tab$Freq > 1,]
  
  duplicates <- degs %>%
    filter(Name %in% dup_tab$Var1) %>%
    group_by(Name) %>%
    summarise(N = n_distinct(Project),
              Project = paste(unique(Project), collapse = ", "),
              Direction = paste(unique(Direction), collapse = ", "),
              .groups = "drop")
  
  # Congruence
  duplicates$Congruence <- sapply(duplicates$Direction, check_congruence)
  duplicates$Congruence[which(duplicates$N == 1)] <- NA
  duplicates <- duplicates[order(duplicates$Congruence),]
  
  # Extract logFC & padj for both comparisons
  sel_logFCA <- degs %>%
    filter(Name %in% duplicates$Name, Project == label1) %>%
    select(Name, log2FC, padj)
  colnames(sel_logFCA) <- c("Name", paste0("log2FC_", label1), paste0("padj_", label1))
  
  sel_logFCB <- degs %>%
    filter(Name %in% duplicates$Name, Project == label2) %>%
    select(Name, log2FC, padj)
  colnames(sel_logFCB) <- c("Name", paste0("log2FC_", label2), paste0("padj_", label2))
  
  sel_logFC <- merge(sel_logFCA, sel_logFCB, by = "Name", all = TRUE)
  res_shared <- merge(duplicates, sel_logFC, by = "Name", all = TRUE)
  res_shared <- res_shared[order(res_shared$N, decreasing = TRUE), ]
  
  df_filt1 <- df_filt1[match(res_shared$Name, df_filt1$Name),]
  df_filt2 <- df_filt2[match(res_shared$Name, df_filt2$Name),]
  
  res_shared[, paste0("log2FC_", label1)] <- df_filt1$log2FC
  res_shared[, paste0("log2FC_", label2)] <- df_filt2$log2FC
  res_shared[, paste0("padj_", label1)] <- df_filt1$pvalue
  res_shared[, paste0("padj_", label2)] <- df_filt2$pvalue
  
  
  # Output folder: W:/mponce/
  comp_degs <- degs[, c(2,15:ncol(degs)-1)] 
  comp_degs <- comp_degs[!duplicated(comp_degs$Name),]
  dim(comp_degs)
  dim(res_shared)
  
  res_counts <- merge(res_shared, comp_degs, by = "Name")
  res_counts$Project <- factor(res_counts$Project, levels = c(paste(label1, label2, sep = ", "), label1, label2), ordered = TRUE)
  res_counts <- res_counts[order(res_counts$Project),]

  # Save data
  addWorksheet(exc, comparison)
  writeData(exc, comparison, res_counts)

  
  
  # VENN DIAGRAM
  # -----------------------
  # All DEGs
  ggvenn(gene_list, 
         digits = 2,                                    # Two decimals
         fill_color = as.vector(color_list$Condition),       # Fill color
         fill_alpha = 0.7,                            # Transparency, default
         stroke_size = 0.25,                          # Line thickness
         set_name_size = 2.5,                           # Default 
         text_size = 2                              # Default
  )
  ggsave(filename = paste("/00_Venn_diagram_", label_plot, ".pdf",  sep = ""), path = dir_fig, height = 4, width = 4, plot = last_plot())
  
  # Upregulated genes
  ggvenn(up_genes, 
         digits = 2,                                    # Two decimals
         fill_color = as.vector(color_list$Condition),       # Fill color
         fill_alpha = 0.7,                            # Transparency, default
         stroke_size = 0.25,                          # Line thickness
         set_name_size = 2.5,                           # Default 
         text_size = 2,                              # Default
  )
  ggsave(filename = paste("/00_Venn_diagram_UP_", label_plot, ".pdf",  sep = ""), path = dir_fig, height = 4, width = 4, plot = last_plot())
  
  # Downregulated genes
  ggvenn(down_genes, 
         digits = 2,                                    # Two decimals
         fill_color = as.vector(color_list$Condition),       # Fill color
         fill_alpha = 0.7,                            # Transparency, default
         stroke_size = 0.25,                          # Line thickness
         set_name_size = 2.5,                           # Default 
         text_size = 2                              # Default
  )
  ggsave(filename = paste("/00_Venn_diagram_DOWN_", label_plot, ".pdf", sep = ""), path = dir_fig, height = 4, width = 4, plot = last_plot())
  
  

  
  # CORRELATION PLOTS
  # ------------------------------------------------------------------------
  
  # Merge data frames
  data_cor <- merge(filtered[which(filtered$Project == label1), c("Comparison", "Name", "DEG", "Direction", "padj", "log2FC")], 
                    filtered[which(filtered$Project == label2), c("Comparison", "Name", "DEG", "Direction", "padj", "log2FC")], 
                    by = "Name", all = TRUE)
  
  # Color label
  data_cor <- data_cor %>% mutate(Colorlab = case_when(
    DEG.x == "YES" & DEG.y == "YES"  ~ "DEG Shared",
    DEG.x == "YES" & DEG.y == "NO" ~ paste("DGEs_", gsub(" ", "_", label1), sep = ""),
    DEG.x == "NO" & DEG.y == "YES"  ~ paste("DGEs_", gsub(" ", "_", label2), sep = ""),
    DEG.x == "NO" & DEG.y == "NO" ~ "Non significant")) 
  
  
  # Correlation plot without the significant events
  data_transc_non <- data_cor[which(data_cor$Colorlab != "Non significant"), ]
  
  # Correlation 
  cor_all <- cor.test(x = data_transc_non$log2FC.x, data_transc_non$log2FC.y, method = "spearman")
  
  ggplot(data_transc_non, aes(x = log2FC.x, y = log2FC.y, color = Colorlab)) +
    geom_point(size = 1.5, alpha = 0.7)+
    scale_color_manual(values = color_list$Correlation)+
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
    labs(x = paste("log2FC", label1), y = paste("log2FC", label2), 
         title = paste(comparison, ": ", label1, " - ", label2, sep = ""),
         subtitle = paste("rho=", round(cor_all$estimate, 2) , "; p=", sprintf("%.2e", cor_all$p.value), sep = ""))+  
    theme(text = element_text(size = 6), axis.title = element_text(size = rel(1.25)))
  ggsave(filename = paste("01_Correlation_", label_plot, ".pdf", sep = ""), height = 5, width = 6, path = dir_fig)
  
  
  # Congruence in the shared genes
  sh_data <- data_cor[which(data_cor$Colorlab == "DEG Shared"), ]
  sh_data <- sh_data %>% mutate(Congruence = case_when(
    Direction.x == "Upregulated" & Direction.y == "Upregulated"  ~ "YES",
    Direction.x == "Upregulated" & Direction.y == "Downregulated"  ~ "NO",
    Direction.x == "Downregulated" & Direction.y == "Upregulated"  ~ "NO",
    Direction.x == "Downregulated" & Direction.y == "Downregulated"  ~ "YES")) 
  
  
  print(table(duplicates$Congruence) == table(sh_data$Congruence))
  
  
  if(nrow(sh_data) >= 3){
    
    # Correlation 
    cor_all <- cor.test(x = sh_data$log2FC.x, sh_data$log2FC.y, method = "spearman")
    
    # Correlation: Check congruence in the shared genes
    ggplot(sh_data, aes(x = log2FC.x, y = log2FC.y, color = Congruence)) +
      geom_point(size = 1.5, alpha = 0.7)+
      scale_color_manual(values =  c("NO" = "red", "YES" = "darkgreen"))+
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
      geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
      labs(x = paste("log2FC", label1), y = paste("log2FC", label2), 
           title = paste(comparison, ": ", label1, " - ", label2, sep = ""),
           subtitle = paste("rho=", round(cor_all$estimate,2) , "; p=", sprintf("%.2e", cor_all$p.value), sep = ""))+  
      theme(text = element_text(size = 6), axis.title = element_text(size = rel(1.25)))
    ggsave(filename = paste("02_Correlation_", label_plot, ".pdf", sep = ""), height = 5, width = 6, path = dir_fig)
    
    # Correlation plot + labels
    ggplot(sh_data, aes(x = log2FC.x, y = log2FC.y, color = Congruence, label = Name)) +
      geom_point(size = 1.5, alpha = 0.7)+
      geom_text_repel(hjust = "outward", max.overlaps = 10, size = 2, na.rm = TRUE, min.segment.length = 0.6)+
      scale_color_manual(values =  c("NO" = "red", "YES" = "darkgreen"))+
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
      geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), col = "gray", linetype = 'dashed')+
      labs(x = paste("log2FC", label1), y = paste("log2FC", label2), 
           title = paste(comparison, ": ", label1, " - ", label2, sep = ""),
           subtitle = paste("rho=", round(cor_all$estimate,2) , "; p=", sprintf("%.2e", cor_all$p.value), sep = ""))+  
      theme(text = element_text(size = 6), axis.title = element_text(size = rel(1.25)))
    ggsave(filename = paste("03_Correlation_", label_plot, ".pdf", sep = ""), height = 5, width = 6, path = dir_fig)
  }
  
  
  
}



############################################################################l####
#                             SAVE DATA
################################################################################

saveWorkbook(exc, file = paste(dir_res, "/Shared_", analysis, "_", project, "_", labs, ".xlsx", sep = ""), overwrite = TRUE)










