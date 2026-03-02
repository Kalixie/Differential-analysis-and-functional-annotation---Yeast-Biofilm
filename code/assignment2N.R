# Load libraries for RNA-seq and visualization ----

library(tximport)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(pheatmap)
library(enrichplot)
library(DOSE)
library(tibble)
library(dplyr)
library(GenomicFeatures)
library(AnnotationDbi)
library(gridExtra)

# Data Preparation ----

# Create sample metadata table

samples <- tribble(
  ~sample, ~stage, ~run,
  "IL20", "Early", "SRR10551665",
  "IL21", "Early", "SRR10551664",
  "IL22", "Early", "SRR10551663",
  "IL23", "Thin", "SRR10551662",
  "IL24", "Thin", "SRR10551661",
  "IL25", "Thin", "SRR10551660",
  "IL29", "Mature", "SRR10551659",
  "IL30", "Mature", "SRR10551658",
  "IL31", "Mature", "SRR10551657"
)

# Set factor levels to define order

samples$stage <- factor(samples$stage,
  levels = c("Early", "Thin", "Mature")
)

# File paths to Salmon quantification outputs

files <- file.path(paste0(samples$run, "_quant"), "quant.sf")

# Name files using sample IDs

names(files) <- samples$sample

# Build transcript database from GTF annotation

txdb <- makeTxDbFromGFF("GCF_000146045.2_R64_genomic.gtf.gz")

k <- keys(txdb, keytype = "TXNAME")

# Create transcript2gene mapping table

tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Import Salmon quantifications

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)


# Run DEQseq2 ----

# Create DESeq dataset object using stage as design factor

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~stage
)

# Run DESeq2 and create plots

dds <- DESeq(dds)

# Comparison pairwise results for each group

res_Mature_vs_Early <- results(dds, contrast = c("stage", "Mature", "Early"))

res_Thin_vs_Early <- results(dds, contrast = c("stage", "Thin", "Early"))

res_Mature_vs_Thin <- results(dds, contrast = c("stage", "Mature", "Thin"))

# Results to data frame and label significant

processresults <- function(res) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)

  df$significant <- ifelse(
    df$padj < 0.05 & abs(df$log2FoldChange) > 1,
    ifelse(df$log2FoldChange > 0, "Up", "Down"),
    "Not Sig"
  )

  na.omit(df)
}

# Process all comparisons

df_MvE <- processresults(res_Mature_vs_Early)

df_TvE <- processresults(res_Thin_vs_Early)

df_MvT <- processresults(res_Mature_vs_Thin)


# Count significantly upregulated genes

countup <- function(df) {
  sum(df$padj < 0.05 & df$log2FoldChange > 1, na.rm = TRUE)
}

# Count significantly downregulated genes

countdown <- function(df) {
  sum(df$padj < 0.05 & df$log2FoldChange < -1, na.rm = TRUE)
}

# Calculate DE gene counts for each comparison

countup(df_MvE)
countup(df_TvE)
countup(df_MvT)

countdown(df_MvE)
countdown(df_TvE)
countdown(df_MvT)

# MA plots for each comparison

plotMA(res_Mature_vs_Early, ylim = c(-5, 5))

plotMA(res_Thin_vs_Early, ylim = c(-5, 5))

plotMA(res_Mature_vs_Thin, ylim = c(-5, 5))

# Extract transformed and normalized counts

vsd <- vst(dds)

# PCA plot for all groups by stage

pcadata <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))

# Plot PCA visualization

ggplot(pcadata, aes(x = PC1, y = PC2, color = stage)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") +
  coord_fixed()

# Volcano plots

plotvolcano <- function(df, title) {
  ggplot(df, aes(log2FoldChange, -log10(pvalue), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 p-value")
}

# Plots for each group pair

p1 <- plotvolcano(df_MvE, "Volcano plot: Mature vs Early")
p2 <- plotvolcano(df_TvE, "Volcano plot: Thin vs Early")
p3 <- plotvolcano(df_MvT, "Volcano plot: Mature vs Thin")

grid.arrange(p1, p2, p3, ncol = 1)

# Function to plot heatmap of top DE genes

plotheatmap <- function(res_object, vsd, dds, n_genes = 20, title = "Top DE Genes") {
  res_clean <- na.omit(res_object)

  top_indices <- head(order(res_clean$padj), n_genes, decreasing = FALSE)
  gene_names <- rownames(res_clean)[top_indices]

  mat <- assay(vsd)[gene_names, ]

  annotation_df <- data.frame(Stage = colData(dds)$stage)
  rownames(annotation_df) <- colnames(mat)

  pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_df,
    main = title,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
}

# Heatmap plot generation

plotheatmap(res_Mature_vs_Early, vsd, dds, title = "Mature vs Early")

plotheatmap(res_Thin_vs_Early, vsd, dds, title = "Thin vs Early")

plotheatmap(res_Mature_vs_Thin, vsd, dds, title = "Mature vs Thin")

# Top DEG tables

getTopDEGs <- function(df, n = 15) {
  df %>%
    arrange(padj) %>%
    mutate(Regulation = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
    dplyr::select(gene, log2FoldChange, padj, Regulation) %>%
    head(n)
}

# Apply to each pairwise

topMvE <- getTopDEGs(df_MvE)
topTvE <- getTopDEGs(df_TvE)
topMvT <- getTopDEGs(df_MvT)

# Run GO ORA and create ----

# Function to run GO enrichment for significant DE genes

runGO <- function(df) {
  sig_genes <- df %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  bg_genes <- df$gene %>%
    na.omit() %>%
    unique()

  enrichGO(
    gene = sig_genes,
    universe = bg_genes,
    OrgDb = org.Sc.sgd.db,
    keyType = "ORF",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
}

GOMvE <- runGO(df_MvE)

GOTvE <- runGO(df_TvE)

GOMvT <- runGO(df_MvT)

# GO analysis plots

dotplot(GOMvE, title = "GO: Mature vs Early")

dotplot(GOTvE, title = "GO: Thin vs Early")

dotplot(GOMvT, title = "GO: Mature vs Thin")


# Upregulated genes

GOup <- function(df) {
  bg_genes <- df$gene %>%
    na.omit() %>%
    unique()
  up_genes <- df %>%
    filter(padj < 0.05 & log2FoldChange > 1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  ego_up <- enrichGO(
    gene = up_genes,
    universe = bg_genes,
    OrgDb = org.Sc.sgd.db,
    keyType = "ORF",
    ont = "BP",
    pvalueCutoff = 0.05
  )
  return(ego_up)
}

# Downregulated genes

GOdown <- function(df) {
  bg_genes <- df$gene %>%
    na.omit() %>%
    unique()
  down_genes <- df %>%
    filter(padj < 0.05 & log2FoldChange < -1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  ego_down <- enrichGO(
    gene = down_genes,
    universe = bg_genes,
    OrgDb = org.Sc.sgd.db,
    keyType = "ORF",
    ont = "BP",
    pvalueCutoff = 0.05
  )
  return(ego_down)
}

GOcompare <- function(df) {
  bg_genes <- df$gene %>%
    na.omit() %>%
    unique()
  up_genes <- df %>%
    filter(padj < 0.05 & log2FoldChange > 1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  down_genes <- df %>%
    filter(padj < 0.05 & log2FoldChange < -1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  compare_df <- compareCluster(
    geneCluster = list(
      Upregulated = up_genes,
      Downregulated = down_genes
    ),
    fun = "enrichGO",
    OrgDb = org.Sc.sgd.db,
    keyType = "ORF",
    ont = "BP"
  )
  return(compare_df)
}

# Results

MvE_up <- GOup(df_MvE)
MvE_down <- GOdown(df_MvE)
MvE_compare <- GOcompare(df_MvE)

dotplot(MvE_up, showCategory = 10, title = "GO: Mature vs Early: Upregulated Genes")

dotplot(MvE_down, showCategory = 10, title = "GO: Mature vs Early: Downregulated Genes")

dotplot(MvE_compare, showCategory = 5, title = "GO: Mature vs Early (Up vs Down regulated)")

TvE_up <- GOup(df_TvE)
TvE_down <- GOdown(df_TvE)
TvE_compare <- GOcompare(df_TvE)

dotplot(TvE_up, showCategory = 10, title = "GO: Thin vs Early: Upregulated Genes")

dotplot(TvE_down, showCategory = 10, title = "GO: Thin vs Early: Downregulated Genes")

dotplot(TvE_compare, showCategory = 5, title = "GO: Thin vs Early (Up vs Down regulated)")

MvT_up <- GOup(df_MvT)
MvT_down <- GOdown(df_MvT)
MvT_compare <- GOcompare(df_MvT)

dotplot(MvT_up, showCategory = 10, title = "GO: Mature vs Thin: Upregulated Genes")

dotplot(MvT_down, showCategory = 10, title = "GO: Mature vs Thin: Downregulated Genes")

dotplot(MvT_compare, showCategory = 5, title = "GO: Mature vs Thin (Up vs Down regulated)")

# Run KEGG ORA and create plots ----

# Run KEGG pathway enrichment on significant genes

keggenrichment <- function(df) {
  sig_genes <- df %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  kegg_enrich <- enrichKEGG(
    gene = sig_genes,
    organism = "sce",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  return(kegg_enrich)
}

# Run KEGG for each comparison

kegg_MvE <- keggenrichment(df_MvE)
kegg_TvE <- keggenrichment(df_TvE)
kegg_MvT <- keggenrichment(df_MvT)

# Visualize KEGG Enrichment Results

dotplot(kegg_MvE, showCategory = 10, title = "KEGG: Mature vs Early")

dotplot(kegg_TvE, showCategory = 10, title = "KEGG: Thin vs Early")

dotplot(kegg_MvT, showCategory = 10, title = "KEGG: Mature vs Thin")

# KEGG compare plots

keggcompare <- function(df) {
  bg_genes <- df$gene %>%
    na.omit() %>%
    unique()

  up_genes <- df %>%
    filter(padj < 0.05 & log2FoldChange > 1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  down_genes <- df %>%
    filter(padj < 0.05 & log2FoldChange < -1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  compare_df <- compareCluster(
    geneCluster = list(
      Upregulated = up_genes,
      Downregulated = down_genes
    ),
    fun = "enrichKEGG",
    organism = "sce",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  return(compare_df)
}

# Comparison of up and downregulated genes

kegg_compare_MvE <- keggcompare(df_MvE)

kegg_compare_TvE <- keggcompare(df_TvE)

kegg_compare_MvT <- keggcompare(df_MvT)

# Plot the KEGG enrichment result comparisons

dotplot(kegg_compare_MvE, showCategory = 5, title = "KEGG: Mature vs Early (Up vs Downregulated)")

dotplot(kegg_compare_TvE, showCategory = 5, title = "KEGG: Thin vs Early (Up vs Downregulated)")

dotplot(kegg_compare_MvT, showCategory = 5, title = "KEGG: Mature vs Thin (Up vs Downregulated)")
