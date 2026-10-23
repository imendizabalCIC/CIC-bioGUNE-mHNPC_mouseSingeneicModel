
#################################################################################
#                           CREATE GENE LENGTH FILE
#################################################################################

# Summary
#---------




# Install libraries
# install.packages("BiocManager")
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "txdbmaker"))

# Load libraries
library(GenomicFeatures)
library(AnnotationDbi)
library(txdbmaker)


# Path Genomes Rocky
path <- "path"


# Create directory
dir.create(file.path(path, "GeneLength"), showWarnings = FALSE)
dir_genelength <- paste(path, "/GeneLength", sep = "")
setwd(dir_genelength)


# Create the gene length files fro human and mouse
for (i in 1:2){
  
  # Load GTF files
  if(i == 1){
    gtf <- paste(path, "/Annotation/Homo_sapiens.GRCh38.110.gtf", sep = "")
    filename <- "Homo_sapiens.GRCh38.110"
  }else{
    gtf <- paste(path, "/Annotation/Mus_musculus.GRCm39.110.gtf", sep = "")
    filename <- "Mus_musculus.GRCm39.110"
  }
  
  # Create a TxDb object from the GTF file
  txdb <- makeTxDbFromGFF(gtf, format = "gtf")
  
  # Load genes data
  txgenes <- genes(txdb, single.strand.genes.only = FALSE)
  
  # Extract gene lengths by summing the lengths of exons for each gene
  gene_length <- sum(width(txgenes))
  
  # Convert to a data frame for easier handling
  gene_length <- data.frame(gene_id = names(gene_length), length = gene_length)
  
  
  
  # Save the gene length file
  write.table(x = gene_length, file = paste("GeneLength_", filename, ".txt", sep = ""), row.names = F, sep = "\t", quote = F)
  
}
