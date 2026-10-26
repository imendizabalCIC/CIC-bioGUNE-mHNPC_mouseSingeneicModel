
# CIC-bioGUNE-mHNPC_mouseSingeneicModel

Code for Cid-Centelles et al., *Cell Death Discovery*, 2026.

## Gene_Expression file

Scripts used for RNA-seq data processing and differential gene expression analysis. Raw sequencing reads were trimmed with Cutadapt[^1], checked for quality using FastQC[^2] and MultiQC[^3], and aligned to the reference genome with STAR[^4]. Differential gene expression analysis was carried out using DESeq2[^5] to identify genes that were differentially expressed between conditions. The folder also contains scripts for quality control, batch effect correction with ComBat-seq[^6], and functional enrichment analysis using GSEA[^7] and g:Profiler[^8].

## mHNPC_WESmouse file

Scripts for the identification of somatic mutations from whole-exome sequencing (WES) data. Raw sequencing reads were aligned to the mouse reference genome GRCm38 using BWA[^9], followed by duplicate marking with Picard[^10] and base quality score recalibration (BQSR) using GATK[^11].

### Single Nucleotide Variants (SNVs)

Somatic SNVs were identified using Mutect2[^12] from GATK with matched normal samples, a mouse germline variant resource, and an interval list corresponding to the WES capture design. Variants were filtered using FilterMutectCalls and annotated with the Ensembl Variant Effect Predictor (VEP)[^13].

### Copy Number Alterations (CNAs)

CNAs were identified using Control-FREEC[^14] with tumor and matched normal BAM files. The analysis was restricted to the targeted exome regions and incorporated GC-content normalization and genome mappability information to reduce technical biases.


[^1]: Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1). https://doi.org/10.14806/ej.17.1.200

[^2]: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[^3]: Ewels P. et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354

[^4]: Dobin A. et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15-21. https://doi.org/10.1093/bioinformatics/bts635

[^5]: Love M.I., Huber W. and Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

[^6]: Zhang Y., Parmigiani G. and Johnson W.E. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*, 2(3), lqaa078. https://doi.org/10.1093/nargab/lqaa078

[^7]: Subramanian A. et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 102(43), 15545-15550. https://doi.org/10.1073/pnas.0506580102

[^8]: Kolberg L. et al. (2023). g:Profiler-interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). *Nucleic Acids Research*, 51(W1), W207-W212. https://doi.org/10.1093/nar/gkad347

[^9]: Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, 25(14), 1754-1760. https://doi.org/10.1093/bioinformatics/btp324

[^10]: Broad Institute (2019). Picard Toolkit. https://broadinstitute.github.io/picard/

[^11]: McKenna A. et al. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*, 20(9), 1297-1303. https://doi.org/10.1101/gr.107524.110

[^12]: Benjamin D. et al. (2019). Calling Somatic SNVs and Indels with Mutect2. *bioRxiv*. https://doi.org/10.1101/861054

[^13]: McLaren W. et al. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17, 122. https://doi.org/10.1186/s13059-016-0974-4

[^14]: Boeva V. et al. (2012). Control-FREEC: a tool for assessing copy number and allelic content using next-generation sequencing data. *Bioinformatics*, 28(3), 423-425. https://doi.org/10.1093/bioinformatics/btr670

