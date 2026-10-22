# CIC-bioGUNE-mHNPC_mouseSingeneicModel

Code for Cid-Centelles et al, Cell Death Discovery, 2006

## mHNPC_WESmouse file: Scripts for the identification of somatic mutations from whole-exome sequencing (WES) data. Raw sequencing reads were aligned to the mouse reference genome GRCm38 using BWA, followed by duplicate marking with **Picard** and base quality score recalibration (BQSR) using GATK.

### Single Nucleotide Variants (SNVs): Somatic SNVs were identified using Mutect2 from GATK with matched normal samples, a mouse germline variant resource, and an interval list corresponding to the WES capture design. Variants were filtered using FilterMutectCalls and annotated with the Ensembl Variant Effect Predictor (VEP).

### Copy Number Alterations (CNAs): CNAs were identified using Control-FREEC with tumor and matched normal BAM files. The analysis was restricted to the targeted exome regions and incorporated GC-content normalization and genome mappability information to reduce technical biases.
