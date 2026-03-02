# Differential analysis and functional annotation - Yeast Biofilm

## Introduction

Through the use of Saccharomyces cerevisiae velum RNA samples, this study aimed to explore the transcriptional changes associated with three different sampling points. Biofilm formation is a survival strategy used by many microorganisms including Saccharomyces cerevisiae, and allows cells to survive harsh environments through altered gene expression (Mardanov et al., 2020). By studying biofilm at the early (38 days), thin formation (83 days), and mature stages (109 days), it is possible to analyze how gene expression changes over time, and what biological processes are associated with said genes. 

The data used in this study were collected from publicly available bulk RNA sequencing data from three stages of velum development (early, thin, and mature) under winemaking conditions, each with three replicate samples (Mardanov et al., 2020). Several computational tools were selected to complete an efficient analysis of the data. 

The pseudoalignment tool Salmon was used for transcript quantification because it provides quick and memory-efficient abundance estimation for RNA-seq data (Patro et al., 2017). Although Salmon cannot provide splice junction information, this was not required for this study, and only gene counts were required for further analysis, so a traditional aligner such as HISAT2 or STAR was not chosen. 

Differential expression analysis was performed in RStudio through the use of the DESeq2 package. DESeq2 uses shrinkage estimation and fold changes to create reliable statistical results, even with the use of smaller sample sizes (Love et al., 2014). DESeq2 was chosen for its widespread use and strong statistical framework. Alternative methods such as EdgeR or limma-voom can also be used for RNA-seq analysis; however would benefit from a larger sample size than the one used for this study (Rapaport et al., 2013).  

To interpret the biological significance of discovered differentially expressed genes, over representation analysis (ORA) using Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) was conducted. ORA identifies categories related to biological processes that occur more frequently through the differentially expressed genes then what is typically expected. Through the usage of ORA and KEGG, functional processes and pathways could be identified as either over or under represented throughout sample comparisons (Young et al., 2010). Overall, the selected workflow provides a reliable and accessible framework for identifying transcriptional and functional changes associated with yeast velum development.

## Methods

### Data acquisition and preprocessing

Raw RNA-seq data for nine samples representing three velum developmental stages (early, thin, and mature) were obtained from a previous study (https://doi.org/10.3389/fmicb.2020.00538). Sequence reads were converted to FASTQ format using SRA tools (v. 3.2.1) (NCBI, 2026). Quality assessment was performed using FastQC (v. 0.12.1) (Andrews, 2010) to evaluate read quality. Although read trimming was performed using fastp (v. 1.1.0) (Chen et al., 2018), post-trimming quality metrics showed no improvement in data therefore, original reads were used for downstream analysis.

### Transcript quantification

Transcript abundance was quantified using Salmon (v. 1.10.3) (Patro et al., 2017). The reference transcriptome for Saccharomyces cerevisiae (R64 assembly) was downloaded from NCBI and indexed using Salmon’s indexing function. Each FASTQ file was then processed with the automatic library flag and the threads flag. The tximport (v. 1.36.1) (Soneson et al., 2015) package in R was used for isolating gene-level counts, and a transcript-to-gene mapping table was generated from the associated GTF file retrieved from NCBI. 

### Differential expression analysis

Gene counts were then analyzed with DESeq2 (v. 1.48.1) (Love et al., 2014). Pairwise comparisons were performed between developmental stages to identify genes significantly up or downregulated during biofilm progression at each stage. Genes were considered differentially expressed if they met thresholds of adjusted p-value < 0.05 and absolute log2 fold change > 1. 

### Data visualization

Principal component analysis (PCA) was performed on count data to examine overall sample clustering and detect outliers. MA plots were generated to showcase global expression changes and identify significantly differentially expressed genes. Volcano plots were used to display the gene expression relationship between each of the groups. Additionally, heatmaps of the top differentially expressed genes were constructed to showcase expression patterns across developmental stages.

### Functional enrichment analysis

Biological interpretation of differentially expressed genes was performed using over-representation analysis (ORA) with the clusterProfiler (v. 4.16.0) (Yu et al., 2012) package. Significant genes were compared against a background set consisting of all genes tested in the differential expression analysis. Gene Ontology and KEGG were used identify biological pathways and processes enriched among upregulated and downregulated genes. Results were visualized using bar plots and dot plots.

## Results

<div align="center">

![pca](https://github.com/Kalixie/Differential-analysis-and-functional-annotation---Yeast-Biofilm/blob/main/figures/pca.png)

</div>

##### Figure #1: PCA Results for Early, Thin, and Mature stages. Clusters can be distinctly identified with PC one accounting for 67% of the variance seen, and PC two accounting for 25%

The results of the Principal component analysis revealed differences in gene expression among the three developmental stages. The first principal component explained 67% of the total variance and clearly separated samples, with Early samples clustering on the negative side and Mature samples clustering on the positive side, with clear distinction (Figure 1). Thin stage samples were positioned between Early and Mature stages, which aligns with it acting as an intermediate developmental state. The second principal component, which accounted for 26% of the variance, further distinguished the Thin samples from the other stages. Replicates within each stage clustered tightly together, indicating groups had low variability within themselves, although the thin group was slightly more spread out. Overall, the PCA demonstrates a difference across stage data that progresses from Early to Thin to Mature.

##### Table #1: Summary of differentially expressed genes across biofilm developmental stages. The table shows the number of genes significantly upregulated, downregulated, and the total number of significant genes for each pairwise comparison.

<div align="center">
  
| Comparison      | Upregulated | Downregulated | Total Significant |
|-----------------|-------------|---------------|-------------------|
| Mature vs Early | 1010        | 862           | 1872              |
| Thin vs Early   | 594         | 532           | 1126              |
| Mature vs Thin  | 731         | 588           | 1319              |

</div>

<div align="center">

![volcano](https://github.com/Kalixie/Differential-analysis-and-functional-annotation---Yeast-Biofilm/blob/main/figures/volcanoplots.png)

</div>

##### Figure #2: Volcano Plot Results for Mature vs Early, Thin vs Early, and Mature vs Thin stages. Log2 Fold Change can be seen on the X-axis, measuring the magnitude of gene expression change. The Y-axis depicts the -log10 p-value, demonstrating the significance levels for each gene.

The results of the volcano plots demonstrated that Mature vs Early showed the largest number of significantly differentially expressed genes at 1872, with large fold changes. The spread across both axes depicts high up and downregulation of genes between early and mature stages. The Thin vs Early also depicted a large number of differentially expressed genes (1126) with a lower magnitude than Mature vs Early, yet similar upregulated (594) and downregulated (532) counts. Mature vs Thin had fewer downregulated genes (588) in comparison to upregulated genes (731) for a total of 1319 significant genes. These volcano plots suggest a transcriptional shift during development, with the largest changes occurring between the Early and Mature stages.

</div>

<div align="center">

![heatmap MvE](https://github.com/Kalixie/Differential-analysis-and-functional-annotation---Yeast-Biofilm/blob/main/figures/Heatmap_MvE.png)

</div>

##### Figure #3: Heatmap Results for the top 20 genes based on padj for the Mature vs Early stage. Red indicates higher expression, while blue indicates lower expression patterns. Gene names are listed on the right, and sample IDs are included below.

The results of the heatmaps generated depict how specific gene logfoldchange values vary between stages for the top 20 differentially expressed genes. One group of genes was highly expressed in Early samples and downregulated in Thin and Mature stages which included genes such as YGR088W and YGR087C. Another group showed the opposite trend, with low expression in Early stages and higher expression in Mature stages, and included genes such as YNR073C and YNR072W. 

</div>

<div align="center">

![GO MvE](https://github.com/Kalixie/Differential-analysis-and-functional-annotation---Yeast-Biofilm/blob/main/figures/GOcomparison_MvE.png)

</div>

##### Figure #4: GO Overrepresentation Analysis Results for the top 5 pathways based on padj for the Mature vs Early stage. The dot size shows the gene ratio, and the y-axis lists the top enriched terms. Color indicates adjusted p-value value. Upregulated and downregulated genes are displayed separately on the x-axis with counts.

Gene Ontology (GO) overrepresentation analysis identified enriched biological processes across pairwise stage comparisons. In the Mature vs Early comparison, upregulated genes were enriched for mitochondrial associated processes, including mitochondrial organization, mitochondrial respiratory chain complex assembly, mitochondrial translation, and energy derivation by oxidation of organic compounds. Downregulated genes were enriched for metabolic and transport related processes, including transmembrane transport, lipid metabolism, monocarboxylic acid metabolism, organic acid metabolism, and oxoacid metabolism. For Mature vs Thin, energy reserve and glycogen metabolic processes were observed to be upregulated, while lipid biosynthetic and metabolic processes and transmembrane transport were downregulated. The upregulation of cytoplasmic translation and ribosome biogenesis was showcased by Thin vs Early, with downregulation of nicotinamide and pyridine nucleotide metabolic processes.

</div>

<div align="center">

![KEGG MvE](https://github.com/Kalixie/Differential-analysis-and-functional-annotation---Yeast-Biofilm/blob/main/figures/KEGGcomparison_MvE.png)

</div>

##### Figure #5: KEGG Overrepresentation Analysis Results for the top 5 pathways based on padj for the Mature vs Early stage. The dot size shows the gene ratio, and the y-axis lists the top related pathways. Color indicates adjusted p-value value. Upregulated and downregulated genes are displayed separately on the x-axis with counts.

For the KEGG analysis, in the Mature vs Early comparison, significantly enriched upregulated pathways included ribosome, proteasome, citrate cycle (TCA cycle), oxidative phosphorylation, and biosynthesis of secondary metabolites. Downregulated pathways were enriched for the biosynthesis of secondary metabolites, fatty acid metabolism, and glycolysis/gluconeogenesis, and carbon metabolism. Similar patterns were observed in the Thin vs Early comparison, although fewer terms reached significance relative to the Mature vs Early stages. Mature vs Thin included upregulated genes related to metabolic processes involving Starch and sucrose, beta-Alanine, and Glycerolipids. Downregulated genes involved processes related to fatty acid and carbon metabolism.   

## Discussion

Biofilm formation, along with related metabolic changes, helps yeast survive harsh conditions such as high ethanol levels and low nutrient availability (Alexandre, 2013). Between the Mature and Early biofilm stages, 1010 upregulated and 862 downregulated genes were identified. This indicates that transcriptional change occurs during biofilm maturation. In contrast, the Thin vs Early and Mature vs Thin comparison showed fewer significant changes, suggesting that the Thin stage represents a transitional biofilm phase between both states.

​
The PCA results further supported this interpretation. Samples clustered strongly by stage, confirming that stage acts as a source of variation in the dataset. Mature and thin samples formed the most distinct clusters, while Thin had a more spread out cluster, suggesting greater variation within that sample group, consistent with its presence as an in-between state of development.

​
Heatmap visualization of the top 20 differentially expressed genes further depicted stage-specific expression patterns. The Mature vs Early comparison showcased fold changes with high contrast between early growth and fully established biofilm states.YGR087C (PDC6) was found to be downregulated in mature samples and more highly expressed in early biofilm stages and codes for an isoenzyme of pyruvate decarboxylase. Pyruvate decarboxylase is involved in the alcoholic fermentation for Saccharomyces cerevisiae, with previous studies identifying its higher presence during the early growth of yeast in ethanol medium (Hohmann, 1991). Upregulated genes between the Mature and Early stages included YNR072W or hexose transporter HXT17. Hxt17 is a transporter for sorbitol and mannitol and its upregulation during mature stages suggests a reliance on transport of alternative carbon sources during stressful conditions (Jordan et al., 2016) (Mardanov et al., 2020).


Overrepresentation analysis provided biological pathway and process information for the transcriptional changes seen through the pairwise differential analysis. Gene Ontology results indicated that many downregulated biological processes between Mature and Early samples were related to transmembrane transport, and metabolic processes relating to monocarboxylic acid, lipid, organic acid, and oxoacid pathways. In contrast, several mitochondrion and oxidation related processes were upregulated in Mature stages. Studies have identified that the deletion of genes involved in mitochondrial organisation led to a decrease in biofilm formation (Vandenbosch et al, 2013). This establishes the importance of the upregulation of said genes for the formation of proper mature biofilms.

​
The KEGG pathway enrichment results included upregulated pathways related to the ribosome, proteasome, TCA cycle, and oxidative phosphorylation. Downregulated pathways included those related to the biosynthesis of secondary metabolites, glycolysis, fatty acid metabolism, and carbon metabolism. Previous findings report that biofilm cells undergo transcriptional reprogramming to support structural matrix changes when forming biofilms (Mundodi et al., 2015). Ribosomes support the production of proteins required for biofilm formation, as demonstrated by the high upregulation identified through KEGG analysis. The enrichment of proteasome pathways indicates that protein quality control is also an essential process during this stage. TCA cycle upregulation is consistent with biofilm development to support growth conditions (Gaupp et al., 2010). 

## References

Alexandre H. (2013). Flor yeasts of Saccharomyces cerevisiae--their ecology, genetics and metabolism. International journal of food microbiology, 167(2), 269–275. https://doi.org/10.1016/j.ijfoodmicro.2013.08.021

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Computer software]. Babraham Bioinformatics. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics (Oxford, England), 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Gaupp, R., Schlag, S., Liebeke, M., Lalk, M., & Götz, F. (2010). Advantage of upregulation of succinate dehydrogenase in Staphylococcus aureus biofilms. Journal of bacteriology, 192(9), 2385–2394.

Hohmann S. (1991). Characterization of PDC6, a third structural gene for pyruvate decarboxylase in Saccharomyces cerevisiae. Journal of bacteriology, 173(24), 7963–7969. https://doi.org/10.1128/jb.173.24.7963-7969.1991

Jordan, P., Choe, J. Y., Boles, E., & Oreb, M. (2016). Hxt13, Hxt15, Hxt16 and Hxt17 from Saccharomyces cerevisiae represent a novel type of polyol transporters. Scientific reports, 6, 23502. https://doi.org/10.1038/srep23502

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Mardanov AV, Eldarov MA, Beletsky AV, Tanashchuk TN, Kishkovskaya SA and Ravin NV (2020) Transcriptome Profile of Yeast Strain Used for Biological Wine Aging Revealed Dynamic Changes of Gene Expression in Course of Flor Development. Front. Microbiol. 11:538. doi: 10.3389/fmicb.2020.00538

Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome Profile of Yeast Strain Used for Biological Wine Aging Revealed Dynamic Changes of Gene Expression in Course of Flor Development. Frontiers in microbiology, 11, 538. https://doi.org/10.3389/fmicb.2020.00538

Mundodi, V., Choudhary, S., Smith, A. D., & Kadosh, D. (2025). Ribosome profiling reveals differences in global translational vs transcriptional gene expression changes during early Candida albicans biofilm formation. Microbiology spectrum, 13(3), e0219524. https://doi.org/10.1128/spectrum.02195-24

National Center for Biotechnology Information (NCBI). (2026). SRA Toolkit. GitHub. https://github.com/ncbi/sra-tools

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

Rapaport, F., Khanin, R., Liang, Y., Pirun, M., Krek, A., Zumbo, P., Mason, C. E., Socci, N. D., & Betel, D. (2013). Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data. Genome biology, 14(9), R95. https://doi.org/10.1186/gb-2013-14-9-r95

Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4, 1521. https://doi.org/10.12688/f1000research.7563.2

Vandenbosch, D., De Canck, E., Dhondt, I., Rigole, P., Nelis, H. J., & Coenye, T. (2013). Genomewide screening for genes involved in biofilm formation and miconazole susceptibility in Saccharomyces cerevisiae. FEMS yeast research, 13(8), 720–730. https://doi.org/10.1111/1567-1364.12071

Young, M.D., Wakefield, M.J., Smyth, G.K. et al. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biol 11, R14 (2010). https://doi.org/10.1186/gb-2010-11-2-r14

Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. Omics : a journal of integrative biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118

