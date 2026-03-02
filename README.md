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

Gene counts were then analyzed with DESeq2. Pairwise comparisons were performed between developmental stages to identify genes significantly up or downregulated during biofilm progression at each stage. Genes were considered differentially expressed if they met thresholds of adjusted p-value < 0.05 and absolute log2 fold change > 1. 

### Data visualization

Principal component analysis (PCA) was performed on count data to examine overall sample clustering and detect outliers. MA plots were generated to showcase global expression changes and identify significantly differentially expressed genes. Volcano plots were used to display the gene expression relationship between each of the groups. Additionally, heatmaps of the top differentially expressed genes were constructed to showcase expression patterns across developmental stages.

### Functional enrichment analysis

Biological interpretation of differentially expressed genes was performed using over-representation analysis (ORA) with the clusterProfiler package. Significant genes were compared against a background set consisting of all genes tested in the differential expression analysis. Gene Ontology Biological Process terms were evaluated to identify biological pathways enriched among upregulated and downregulated genes. Results were visualized using bar plots and dot plots.

## Results

The results of the Principal component analysis revealed large differences in gene expression among the three developmental stages. The first principal component explained 67% of the total variance and clearly separated samples, with Early samples clustering on the negative side and Mature samples clustering on the positive side, with great distinction. Thin stage samples were positioned between Early and Mature stages, which aligns with it acting as an intermediate developmental state. The second principal component, which accounted for 26% of the variance, further distinguished the Thin samples from the other stages. Replicates within each stage clustered tightly together, indicating groups had low variability within themselves, although the thin group was slightly more spread out. Overall, the PCA demonstrates a difference across stage data that follows a clear progression from Early to Thin to Mature.

##### Table #1: Summary of differentially expressed genes across biofilm developmental stages. The table shows the number of genes significantly upregulated, downregulated, and the total number of significant genes for each pairwise comparison.
<div align="center">
| Comparison      | Upregulated | Downregulated | Total Significant |
|-----------------|-------------|---------------|-------------------|
| Mature vs Early | 1010        | 862           | 1872              |
| Thin vs Early   | 594         | 532           | 1126              |
| Mature vs Thin  | 731         | 588           | 1319              |
</div>

The results of the volcano plots demonstrated that Mature vs Early showed the largest number of significantly differentially expressed genes at 1872, with large fold changes. The widespread along both axes showcases substantial transcriptional changes between early and mature stages. The Thin vs Early also depicted a large number of differentially expressed genes (1126) with a lower magnitude than Mature vs Early, yet similar upregulated (594) and downregulated (532) counts. Mature vs Thin had fewer downregulated genes (588) in comparison to upregulated genes (731) for a total of 1319 significant genes. These volcano plots suggest a transcriptional shift during development, with the largest changes occurring between the Early and Mature stages.

The results of the heatmaps generated depict how specific gene logfoldchange values vary between stages for the top 20 differentially expressed genes. One group of genes was highly expressed in Early samples and downregulated in Thin and Mature stages which included genes such as YGR088W and YGR087C. Another group showed the opposite trend, with low expression in Early stages and higher expression in Mature stages, and included genes such as YNR073C and YNR072W. 

Gene Ontology (GO) overrepresentation analysis identified enriched biological processes across pairwise stage comparisons. In the Mature vs Early comparison, upregulated genes were enriched for mitochondrial associated processes, including mitochondrial organization, mitochondrial respiratory chain complex assembly, mitochondrial translation, and energy derivation by oxidation of organic compounds. Downregulated genes were enriched for metabolic and transport related processes, including transmembrane transport, lipid metabolism, monocarboxylic acid metabolism, organic acid metabolism, and oxoacid metabolism. For Mature vs Thin, energy reserve and glycogen metabolic processes were observed to be upregulated, while lipid biosynthetic and metabolic processes and transmembrane transport were downregulated. The upregulation of cytoplasmic translation and ribosome biogenesis was showcased by Thin vs Early, with downregulation of nicotinamide and pyridine nucleotide metabolic processes.

For the KEGG analysis, in the Mature vs Early comparison, significantly enriched upregulated pathways included ribosome, proteasome, citrate cycle (TCA cycle), oxidative phosphorylation, and biosynthesis of secondary metabolites. Downregulated pathways were enriched for the biosynthesis of secondary metabolites, fatty acid metabolism, and glycolysis/gluconeogenesis, and carbon metabolism. Similar patterns were observed in the Thin vs Early comparison, although fewer terms reached significance relative to the Mature vs Early contrast. Mature vs Thin included upregulated genes related to metabolic processes such as Starch and sucrose, beta-Alanine, Glycerolipid. Downregulated included fatty acid and carbon metabolism.   

## Discussion



