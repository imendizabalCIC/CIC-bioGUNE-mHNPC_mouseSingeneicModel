################################################################################
#| POST-CONTROL-FREEC CNA GENE OVERLAP ANALYSIS AND VISUALIZATION (MOUSE WES)
################################################################################
#| Date: 02/10/2025
#| Author: Ivana Rondon Lorefice
#|
#| Description:
#|
#| Downstream analysis of copy number alterations (CNAs) detected by CONTROL-FREEC 
#| in mouse whole exome sequencing (WES) data.  
#| - Reads per-sample CONTROL-FREEC CNA calls and filters events by size.  
#| - Maps CNAs to Ensembl gene annotations (mmusculus v102) and retains overlaps with 
#|   protein-coding genes.  
#| - Aggregates CNA–gene overlaps across samples into a unified dataset.  
#| - Produces summary statistics and multiple plots for visualization of CNA 
#|   distribution, oncogene enrichment, and sample-level CNA burden.  
#|
#| Workflow:
#|
#| Setup
#|   - Load libraries (GenomicRanges, maftools, GenVisR, ggplot2, etc.).  
#|   - Set working directory and define sample list from VCF indices.  
#|   - Retrieve Ensembl annotations (gene biotypes, positions) + homologs to human genes.  
#|
#| CNA–gene overlap
#|   - Read CONTROL-FREEC CNA results (per sample).  
#|   - Compute CNA widths; filter out very small (<100 bp) and very large (>300 kb) events.  
#|   - Create GRanges for CNAs and genes; identify overlaps.  
#|   - Retain protein-coding overlaps and combine into a master table.  
#|
#| Visualizations
#|   - Gene-level CNA counts per sample (barplots).  
#|   - CNA frequency plot (per chromosome, cnFreq, mm10).  
#|   - CNA size distributions (density + histogram).  
#|   - Genes per CNA distributions (density + histogram).  
#|   - Oncogene enrichment analysis (mouse ↔ human homologs).  
#|   - CNA burden per sample (fraction of exome affected).  
#|
#| Outputs:
#|   - Results/combined_data.xlsx                          # merged per-gene CNA calls  
#|   - Results/Plots/Barplot_Genes-Samples_<sample>.*      # barplots of affected genes  
#|   - Results/Plots/cnFreq_plot_<sample>.*                # CNA frequency plots  
#|   - Results/Plots/Density/Histogram_*_width_kb_<sample>.* # CNA size distributions  
#|   - Results/Plots/Density/Histogram_*_genes_per_CNA_<sample>.* # genes per CNA  
#|   - Results/Plots/Barplot_Oncogenes-Samples_<sample>.*  # oncogene proportion plots  
#|   - Results/Plots/CNA_Type-Gain_loss_Bar_<sample>.*     # oncogene CNA class breakdown  
#|   - Results/Plots/Barplot_CNA_Burden_<sample>.*         # CNA burden per sample  
#|
#| Comments:
#|   - Tailored for CONTROL-FREEC CNA outputs in mouse WES data (mm10).  
#|   - Filters and cutoffs (CNA width, oncogene list, exonic length) are configurable.  
#|   - Example plots are shown for sample "PL051", but workflow generalizes to all samples.  
################################################################################


################################################################################+
#| LIBRARY AND DATA
################################################################################
packages <- c(
  "ComplexHeatmap", "ggplot2", "GenomicRanges", "gprofiler2",
  "IRanges", "pheatmap", "biomaRt", "TxDb.Mmusculus.UCSC.mm39.knownGene",
  "viridis", "GenVisR", "absCNseq", "dplyr", "maftools",
  "ggsci", "VennDiagram", "RColorBrewer", "CODEX2",
  "BSgenome.Mmusculus.UCSC.mm10", "tidyr", "dotenv"
)

invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}))

#| For plots
theme_set(theme_classic())

#| Setting working directory
working_Dir <- "~/05_CNAs_CONTROL_FREEC_MP/"
setwd(working_Dir)

#| Reading mouse information (containing gene type info)
ensembl_grcm38 <- useEnsembl(biomart = "genes", 
                             dataset = "mmusculus_gene_ensembl", 
                             version = 102)
mouse_annotation_gene_type <- getBM(attributes = c("gene_biotype", "external_gene_name", 
                                                   "chromosome_name", "start_position", "end_position"),
                                    mart = ensembl_grcm38)
names(mouse_annotation_gene_type) <- c("gene_type","gene_name", "CONTIG", "START", "END" )
mouse_annotation_gene_type$CONTIG <- paste("chr",mouse_annotation_gene_type$CONTIG, sep ="")
mouse_annotation_gene_type$CONTIG <- gsub("chr", "", mouse_annotation_gene_type$CONTIG)

#| Reading sample names
Samples <- list.files(path = "X:/irondon/AC-82_WESmouse/3_SNVs_GATK_MP/3_SNVs_step1_TUMOR/", pattern = "*_filtered_sorted.vcf.gz.tbi")
Samples <- gsub("_filtered_sorted.vcf.gz.tbi", "", Samples)
################################################################################


################################################################################
#| GETTING CNAs FROM CONTROL FREEC OUTPUT
################################################################################
#| Minimum width of CNA
min_width <- 300000

#| To save results in a list for every sample
overlap_data_list <- list()

#| Run a loop 
for (i in 1:length(Samples)){
  
  #| Selecting the first sample
  df <- read.table(paste("Results/",Samples[i], "_mapped.rg.sorted.dup.bqsr.bam_CNVs", sep =""), header =F)
  names(df) <- c("CONTIG", "START","END", "COPYNUMBER", "TYPE")
  df$Sample_id <- NA
  df$Sample_id <- Samples[i]
  
  #| Width of the CNA
  df$width_gene <- abs(df$START - df$END)
  
  df <- df[which(df$width_gene > 100),]
  df <- df[which(df$width_gene < min_width),]
  
  #| Create CNA GRanges object
  cna_gr <- GRanges(
    seqnames = df$CONTIG,        
    ranges = IRanges(start = df$START, end =df$END),
    copynumber = df$COPYNUMBER,
    type = df$TYPE,
    width_gene =df$width_gene,
    sample_id = df$Sample_id
  )
  
  #| Create Gene GRanges object
  gene_gr <- GRanges(
    seqnames = mouse_annotation_gene_type$CONTIG,        # Chromosome/Contig information
    ranges = IRanges(start = mouse_annotation_gene_type$START, end = mouse_annotation_gene_type$END),  # Start and end positions
    gene_name = mouse_annotation_gene_type$gene_name,
    gene_type = mouse_annotation_gene_type$gene_type
  )

  
  #| Find overlaps
  overlaps <- findOverlaps(cna_gr, gene_gr)
  
  overlap_data <- data.frame(chr = seqnames(cna_gr[queryHits(overlaps)]),
                             cna_index = queryHits(overlaps),
                             cna_width = width(cna_gr[queryHits(overlaps)]),
                             gene_name = mcols(gene_gr[subjectHits(overlaps)])$gene_name,
                             gene_type = mcols(gene_gr[subjectHits(overlaps)])$gene_type,
                             cna_start = start(cna_gr[queryHits(overlaps)]),
                             cna_end = end(cna_gr[queryHits(overlaps)]),
                             sample_id = mcols(cna_gr[queryHits(overlaps)])$sample_id,
                             copynumber = mcols(cna_gr[queryHits(overlaps)])$copynumber,
                             type = mcols(cna_gr[queryHits(overlaps)])$type
                             )
  
  overlap_data <- overlap_data[which(overlap_data$gene_type == "protein_coding"),]
  overlap_data_list[[Samples[i]]] <- overlap_data

}

#| Combining and saving data
combined_data <- do.call(rbind, overlap_data_list)
combined_data$sample_id[which(combined_data$sample_id == "T051")] <- "PL051"
writexl::write_xlsx(combined_data, "Results/combined_data.xlsx")
combined_data_PL051 <- combined_data[which(combined_data$sample_id =="PL051"),]
################################################################################

################################################################################
#| PLOTS
################################################################################

#####| ------  Barplot of gene number per sample -------------------------------

ggplot(combined_data_PL051, aes(x = sample_id, fill = type)) +
  geom_bar() +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold")
  ) +
  xlab("") +
  ylab("Number of genes") +
  labs(fill ="CNA Class")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))
ggsave("Results/Plots/Barplot_Genes-Samples_PL051.pdf", height = 6, width = 6,units  = "cm")
ggsave("Results/Plots/Barplot_Genes-Samples_PL051.png", height = 6, width = 6,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Barplot_Genes-Samples_PL051.svg", width = 6, height = 6, units = "cm")


#####| ------ CNA Frequency ----------------------------------------------------

#| PL051
new_df <- combined_data[,c("chr","cna_start", "cna_end", "copynumber", "sample_id")]
names(new_df) <- c("chromosome", "start", "end", "segmean", "sample")
new_df <- new_df[which(new_df$sample == "PL051"),]
new_df$chromosome <- gsub("chr", "", new_df$chromosome)
new_df$segmean <- as.integer(new_df$segmean)

table(combined_data[,c("chr", "type")])
cnFreq(new_df, genome="mm10",
       CN_Loss_colour="#56B4E9",
       CN_Gain_colour ="#E69F00",
       CN_low_cutoff = 2,
       CN_high_cutoff = 2) +
  theme_classic() +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold")
  )
ggsave("Results/Plots/cnFreq_plot_PL051.pdf", height = 6, width = 14.5,units  = "cm")
ggsave("Results/Plots/cnFreq_plot_PL051.png", height = 6, width = 14.5,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/cnFreq_plot_PL051.svg", height = 6, width = 14.5, units = "cm")


#####| ------ Width per CNA ----------------------------------------------------

combined_data_PL051$cna_width_kb <- combined_data_PL051$cna_width/100
ggplot(combined_data_PL051, aes(x = cna_width_kb, color = sample_id, fill = sample_id)) +
  geom_density(alpha = 0.3) +  # add transparency for fill
  scale_fill_manual(values = c("#00441B")) +
  scale_color_manual(values = c("#00441B")) +
  labs(x = "CNA Width (kb)",
       y = "Density",
       color = "Sample ID",
       fill = "Sample ID") +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"), legend.position = "none"
  )
ggsave("Results/Plots/Density_plot_width_kb_PL051.pdf", height = 6, width = 9,units  = "cm")
ggsave("Results/Plots/Density_plot_width_kb_PL051.png", height = 6, width = 9,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Density_plot_width_kb_PL051.svg",  height = 6, width = 9, units = "cm")

ggplot(combined_data_PL051, aes(x = cna_width_kb, color = sample_id, fill = sample_id)) +
  geom_histogram(alpha = 0.3, binwidth = 2) +  
  scale_fill_manual(values = c("#00441B")) +
  scale_color_manual(values = c("#00441B")) +
  labs(x = "CNA Width (kb)",
       y = "Frequency",
       color = "Sample ID",
       fill = "Sample ID") +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"), legend.position = "none") 
ggsave("Results/Plots/Histogram_plot_width_kb_PL051.pdf", height = 6, width = 9,units  = "cm")
ggsave("Results/Plots/Histogram_plot_width_kb_PL051.png", height = 6, width = 9,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Histogram_plot_width_kb_PL051.svg",  height = 6, width = 9, units = "cm", dpi = 300)



#####| ------ Genes per CNA ----------------------------------------------------
per_cna <- combined_data%>%
  group_by(sample_id, chr, cna_index, cna_start, cna_end, type, copynumber) %>%
  summarise(
    n_genes = n_distinct(gene_name),
    genes   = paste(sort(unique(gene_name)), collapse = ";"),
    .groups = "drop"
  )


ggplot(per_cna[which(per_cna$sample_id =="PL051"),], aes(x = n_genes, color = sample_id, fill = sample_id)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("#00441B")) +
  scale_color_manual(values = c("#00441B")) +
  labs(x = "Number of Genes per CNA",
       y = "Density",
       color = "Sample ID",
       fill = "Sample ID") +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"), legend.position = "none"
  )
ggsave("Results/Plots/Density_plot_number_of_genes_per_CNA_PL051.pdf", height = 6, width = 9,units  = "cm")
ggsave("Results/Plots/Density_plot_number_of_genes_per_CNA_PL051.png", height = 6, width = 9,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Density_plot_number_of_genes_per_CNA_PL051.svg", height = 6, width = 9, units = "cm")


ggplot(per_cna[which(per_cna$sample_id =="PL051"),], aes(x = n_genes, color = sample_id, fill = sample_id)) +
  geom_histogram(alpha = 0.3, binwidth = 1) +
  scale_fill_manual(values = c("#00441B")) +
  scale_color_manual(values = c("#00441B")) +
  labs(x = "Number of Genes per CNA",
       y = "Frequency",
       color = "Sample ID",
       fill = "Sample ID") +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"), legend.position =  "none")
ggsave("Results/Plots/Histogram_plot_number_of_genes_per_CNA_PL051.pdf", height = 6, width = 9,units  = "cm")
ggsave("Results/Plots/Histogram_plot_number_of_genes_per_CNA_PL051.png", height = 6, width = 9,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Histogram_plot_number_of_genes_per_CNA_PL051.svg", height = 6, width = 9, units = "cm")


#####| ------ Obtaining the proportions of oncogenes ---------------------------
oncogenes <- read.csv("Oncogene_Gene_Name_List.csv")

genome<-"mmusculus_gene_ensembl"
version<-"102"
ensembl <- useEnsembl(biomart = "genes", dataset = genome, version = version)
annotations <- getBM(attributes = c("gene_biotype", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart = ensembl)
names(annotations) <- c("gene_type","gene_name", "CONTIG", "START", "END" )
annotations$CONTIG <- paste0("chr",annotations$CONTIG)

homolog_info <- getBM(
  attributes = c("external_gene_name",  # Mouse Ensembl ID (for merging)
                 "hsapiens_homolog_associated_gene_name", 
                 "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_orthology_type"),
  mart = ensembl
)
names(homolog_info) <- c(c("gene_name",  # Mouse Ensembl ID (for merging)
                           "hsapiens_homolog_associated_gene_name", 
                           "hsapiens_homolog_ensembl_gene", 
                           "hsapiens_homolog_orthology_type"))
annotations_homolog <- merge(annotations, homolog_info, by = "gene_name")
annotations_homolog <- annotations_homolog[which(annotations_homolog$hsapiens_homolog_associated_gene_name !=""),]

oncogenes_mouse_gene_name <- annotations_homolog$gene_name[which(annotations_homolog$hsapiens_homolog_associated_gene_name %in%  oncogenes$gene_name)]


#####| ------ Number of Homolog genes ------------------------------------------
length(oncogenes_mouse_gene_name)

PL051_sample <- combined_data[which(combined_data$sample_id == "PL051"),]

df_data <- data.frame(sample = c("PL051"),
                      fraction = c( (length(PL051_sample$gene_name[which(PL051_sample$gene_name %in% oncogenes_mouse_gene_name)]))*100/length(oncogenes_mouse_gene_name)))


ggplot(df_data, aes(x = reorder(sample, fraction), y = fraction, fill = sample)) +
  geom_bar(stat = "identity") +
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Percentage of Oncogenes") +
  ylim(0,10)+
  scale_fill_manual(values = c("#00441B"))
ggsave("Results/Plots/Barplot_Oncogenes-Samples_PL051.pdf", height = 6, width =5,units  = "cm")
ggsave("Results/Plots/Barplot_Oncogenes-Samples_PL051.png", height = 6, width = 5,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Barplot_Oncogenes-Samples_PL051.svg", height = 6,  width = 5,units = "cm")



#####| ------ Loss and Gain per sample -----------------------------------------
# Prepare data
cna_data <- rbind(PL051_sample[which(PL051_sample$gene_name %in% oncogenes_mouse_gene_name),])
cna_data$direction <- ifelse(cna_data$type == "loss", -1, 1)

# Plot PL051
ggplot(cna_data, aes(x = gene_name, y = direction, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  facet_wrap(~sample_id) +
  scale_y_continuous(breaks = c(-1, 1), labels = c("Loss", "Gain")) +
  labs(
       x = "Gene", y = "CNA Class") +
  theme_minimal() +
  scale_fill_manual(values = c( "#E69F00", "#56B4E9"))+
  theme(
    text = element_text(size = 10, family = "sans",face = "bold"), 
    axis.text.x = element_text(color = "black",angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("Results/Plots/CNA_Type-Gain_loss_Bar_PL051.pdf", height = 6, width =6,units  = "cm")
ggsave("Results/Plots/CNA_Type-Gain_loss_Bar_PL051.png", height = 6, width = 6,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/CNA_Type-Gain_loss_Bar_PL051.svg", height = 6, width = 6, units = "cm")


#####| ------ CNA burden -------------------------------------------------------
total_exonic_length <- 37546581
Samples <- c("PL051")
cnv_burden <- c()
suma <- c()

for (sample in Samples){
  
  unique_cna_width <- combined_data[which(combined_data$sample_id == sample),] %>%
    distinct(cna_index, cna_width)
  
  cnv_burden <- c(cnv_burden, sum(unique_cna_width$cna_width)/total_exonic_length)
  suma <- c(suma,sum(unique_cna_width$cna_width) )
}

df <- data.frame(sample = Samples,
                 cnv_burden= cnv_burden)

ggplot(df, aes(x = sample,y =cnv_burden, fill = sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#00441B"))+
  theme(
    text = element_text(size = 10, family = "sans"), 
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(linewidth = 0.176),
    axis.line.y = element_line(linewidth = 0.176),
    axis.ticks  = element_line(linewidth = 0.176),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "bold"),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("CNA Burden") +
  ylim(0,1)
ggsave("Results/Plots/Barplot_CNA_Burden_PL051.pdf", height = 6, width = 5,units  = "cm")
ggsave("Results/Plots/Barplot_CNA_Burden_PL051.png", height = 6, width = 5,units  = "cm",dpi = 300, bg = "white" )
ggsave("Results/Plots/Barplot_CNA_Burden_PL051.svg",height = 6, width = 5,  units = "cm")
################################################################################

