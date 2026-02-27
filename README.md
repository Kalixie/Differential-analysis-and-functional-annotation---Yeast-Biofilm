# Differential analysis and functional annotation - Yeast Biofilm

## Introduction

Through the use of Saccharomyces cerevisiae velum RNA samples, this study aimed to explore the transcriptional changes associated with three different sampling points. Biofilm formation is a survival strategy used by many microorganisms including Saccharomyces cerevisiae, and allows cells to survive harsh environments through altered gene expression (Mardanov et al., 2020). By studying biofilm at the early, thin formation, and mature stages, it is possible to analyze how gene expression changes over time, and what biological processes are associated with said genes. 

The data used in this study were collected from publicly available bulk RNA sequencing data from three stages of velum development (early, thin, and mature), each with three replicates. 
Several computational tools were selected to complete an efficient analysis of the data. 

The pseudoalignment tool Salmon was used for transcript quantification because it provides quick and memory-efficient abundance estimation for RNA-seq data (Patro et al., 2017). Although Salmon cannot provide splice junction information, this was not required for this study, and only gene counts were required for further analysis, so a traditional aligner such as HISAT2 or STAR was not chosen. 

Differential expression analysis was performed in RStudio through the use of the DESeq2 package. DESeq2 uses shrinkage estimation and fold changes to create reliable statistical results, even with the use of smaller sample sizes (Love et al., 2014). DESeq2 was chosen for its widespread use and strong statistical framework. Alternative methods such as EdgeR or limma-voom can also be used for RNA-seq analysis; however would benefit from a larger sample size than the one used for this study (Rapaport et al., 2013).  

To interpret the biological significance of discovered differentially expressed genes, over representation analysis (ORA) using Gene Ontology (GO) terms was conducted. ORA identifies categories related to biological processes that occur more frequently through the differentially expressed genes then what is typically expected. Through the usage of ORA, functional processes could be identified as either over or under represented throughout sample comparisons (Young et al., 2010). Overall, the selected workflow provides a reliable and accessible framework for identifying transcriptional and functional changes associated with yeast velum development.

## Methods

### Data acquisition and preprocessing

Raw RNA-seq data for nine samples representing three velum developmental stages (early, thin, and mature) were obtained from a previous study (https://doi.org/10.3389/fmicb.2020.00538). Sequence reads were converted to FASTQ format using SRA tools (v. 3.2.1). Quality assessment was performed using FastQC (v. 0.12.1) to evaluate read quality. Although read trimming was performed using fastp (v. 1.1.0), post-trimming quality metrics showed no improvement in data therefore, original reads were used for downstream analysis.

### Transcript quantification

Transcript abundance was quantified using Salmon (v. 1.10.3). The reference transcriptome for Saccharomyces cerevisiae (R64 assembly) was downloaded from NCBI and indexed using Salmon’s indexing function. Each FASTQ file was then processed with the automatic library flag and the threads flag. The tximport package in R was used for isolating gene-level counts, and a transcript-to-gene mapping table was generated from the associated GTF file retrieved from NCBI. 

### Differential expression analysis



### Read Acquisition 
