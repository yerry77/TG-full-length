# Description: This script processes DEG results and DET results to generate a GO enrichment network plot.
# DEG : Genes differentially expressed in the control group and allergic rhinitis in TG
# DET : Genes with differentially expressed transcripts(isoforms) in the TG control group and allergic rhinitis

# GO enrichment network plot of DEG results
#  Load required libraries 
library(ggplot2)
library(dplyr)
library(ggrastr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggrepel)
library(patchwork)
library(readr)
library(clusterProfiler)
library(aPEAR)

# GO enrichment network plot of DEG results
#  Step 1: Read DEG results 
#  The DEG results file should contain columns such as gene, log2FoldChange, pvalue, padj
df <- read.table("/path/to/DEG_results.tsv",
                 header = TRUE, sep = "\t")

#  Step 2: Read novel isoform information 
#  This file contains information on genes that have novel isoforms
novel_isoforms_file <- "/path/to/novel_isoforms_info.tsv"
novel_isoforms <- read_tsv(novel_isoforms_file, show_col_types = FALSE)

novel_gene_ids <- unique(novel_isoforms$gene_id)
gene_ids_with_novel_isoforms <- unique(novel_isoforms$gene_id)

# # Add a new column 'gene_type' to indicate whether a DEG has a novel isoform
df <- df %>%
  mutate(gene_type = ifelse(df$gene %in% novel_gene_ids, 
                            "with_novel_isoform", "without_novel_isoform"))

# Clean gene_id (remove version number)
df$gene_id_clean <- sub("\\..*$", "", df$gene)

# Map ENSEMBL IDs to gene symbols
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = df$gene_id_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
df$gene_symbol <- gene_symbols

#  Step 3: Add group, color, log10_pvalue, and gene_type 
#  Group genes by significance and fold change
df <- df %>%
  mutate(
    group = case_when(
      pvalue < 0.05 & log2FoldChange > 0.3  ~ "Upregulated",
      pvalue < 0.05 & log2FoldChange < -0.3 ~ "Downregulated",
      TRUE                                  ~ "Not Significant"
    ),
    color = case_when(
      group == "Upregulated"   ~ "#ED0000FF",
      group == "Downregulated" ~ "#00468BFF",
      TRUE                     ~ "grey"
    ),
    log10_pvalue = -log10(pvalue),
    gene_type = ifelse(df$gene %in% gene_ids_with_novel_isoforms,
                       "with_novel_isoform", "without_novel_isoform")
  )

# Filter invalid padj values
df <- df %>% filter(!is.na(padj) & padj > 0)

# Filter significant genes
df_significant <- df %>% filter(group != "Not Significant")

#  Step 4: Map significant DEG symbols to ENTREZ IDs for GO enrichment
gene <- unique(df_significant$gene_symbol)
symbol2entrez <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

#  Step 5: GO enrichment for significant genes 
ego <- enrichGO(
  gene          = symbol2entrez$ENTREZID, # ENTREZ IDs of significant DEGs
  OrgDb         = org.Mm.eg.db,           # Mouse annotation database
  keyType       = "ENTREZID",             # Type of input gene IDs
  ont           = "ALL",                  # Ontology: BP, CC, MF, or ALL
  pAdjustMethod = "BH",                   # Multiple testing correction
  pvalueCutoff  = 0.05,                   # p-value threshold
  qvalueCutoff  = 0.05,                   # q-value threshold
  readable      = TRUE                     # Convert back to gene SYMBOLs
)

# Filter GO terms with adjusted p-value < 1e-2
x.data.to_plot <- ego@result[ego@result$p.adjust < 1e-2,]
x.data.to_plot.BP <- x.data.to_plot[x.data.to_plot$ONTOLOGY == "BP",]
x.data.to_plot.CC <- x.data.to_plot[x.data.to_plot$ONTOLOGY == "CC",]
x.data.to_plot.MF <- x.data.to_plot[x.data.to_plot$ONTOLOGY == "MF",]

#  Step 6: Plot enrichment networks 
p1 <- enrichmentNetwork(x.data.to_plot.BP, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-5, -2), direction = 1)

p2 <- enrichmentNetwork(x.data.to_plot.CC, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-5, -2), direction = 1)

p3 <- enrichmentNetwork(x.data.to_plot.MF, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-5, -2), direction = 1)

p4 <- enrichmentNetwork(x.data.to_plot, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-5, -2), direction = 1)

#  Step 7: Save plots and results 
output_folder <- "/path/to/output_folder"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

pdf(file.path(output_folder, "DEG_GO_enrichment_network.pdf"), width = 10, height = 10)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
write.csv(ego@result, file.path(output_folder, "DEG_GO_enrichment_network.csv"), row.names = FALSE)

# GO enrichment for significant genes with DET 
#  Step 8: DET isoform analysis 
#  The DET file should contain isoform IDs, log2FoldChange, pvalue, padj
df_isoform <- read.table("/path/to/DET_isoforms_results.tsv",
                         header = TRUE, sep = "\t")

# Read novel isoform info
novel_isoforms_file <- "/path/to/novel_isoforms_info.tsv"
novel_isoforms <- read_tsv(novel_isoforms_file, show_col_types = FALSE)
novel_isoforms_data <- novel_isoforms$isoform

# Filter invalid p-values
df_isoform <- df_isoform %>% filter(!is.na(pvalue) & pvalue > 0)

# Group by significance and fold change, define colors, annotate novelty
df_isoform <- df_isoform %>%
  mutate(
    group = case_when(
      pvalue < 0.05 & log2FoldChange > 0.3  ~ "Upregulated",
      pvalue < 0.05 & log2FoldChange < -0.3 ~ "Downregulated",
      TRUE                                  ~ "Not Significant"
    ),
    color = case_when(
      group == "Upregulated"   ~ "#ED0000FF",
      group == "Downregulated" ~ "#00468BFF",
      TRUE                     ~ "grey"
    ),
    log10_pvalue = -log10(pvalue),
    is_novel = ifelse(isoform %in% novel_isoforms_data, "Novel", "Known")
  )

# Join with symbol gene info
# Read additional isoform annotation file
symbol_gene_info <- "/path/to/all_isoforms_info.tsv"
symbolgene_info <- read_tsv(symbol_gene_info, show_col_types = FALSE)
df_isoform <- df_isoform %>% left_join(symbolgene_info, by = "isoform")

# Filter significant isoforms
df_isoform_significant <- df_isoform %>% filter(group != "Not Significant")
df_isoform_significant_gene <- unique(df_isoform_significant$symbol_gene)

# Map symbols to ENTREZ IDs for enrichment
symbol2entrez2 <- bitr(df_isoform_significant_gene, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO enrichment for DET genes with novel isoforms
ego3 <- enrichGO(
  gene          = symbol2entrez2$ENTREZID, # ENTREZ IDs of significant isoform genes
  OrgDb         = org.Mm.eg.db,            # Mouse annotation database
  keyType       = "ENTREZID",              # Type of input gene IDs
  ont           = "ALL",                   # Ontology: BP, CC, MF, or ALL
  pAdjustMethod = "BH",                    # Multiple testing correction
  pvalueCutoff  = 0.05,                    # p-value threshold
  qvalueCutoff  = 0.05,                    # q-value threshold
  readable      = TRUE                      # Convert back to gene SYMBOLs
)

# Filter by adjusted p-value < 1e-2
x.data.to_plot3 <- ego3@result[ego3@result$p.adjust < 1e-2,]
x.data.to_plot3.BP <- x.data.to_plot3[x.data.to_plot3$ONTOLOGY == "BP",]
x.data.to_plot3.CC <- x.data.to_plot3[x.data.to_plot3$ONTOLOGY == "CC",]
x.data.to_plot3.MF <- x.data.to_plot3[x.data.to_plot3$ONTOLOGY == "MF",]

# Plot enrichment networks
p5 <- enrichmentNetwork(x.data.to_plot3.BP, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-10, -2), direction = 1)

p6 <- enrichmentNetwork(x.data.to_plot3.CC, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-10, -2), direction = 1)

p7 <- enrichmentNetwork(x.data.to_plot3.MF, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-10, -2), direction = 1)

p8 <- enrichmentNetwork(x.data.to_plot3, drawEllipses = TRUE, fontSize = 2.5,
                        colorBy = "p.adjust", colorType = "pval", nodeSize="Count",
                        minClusterSize = 1, pCutoff = -20) +
      scale_color_viridis(limits = c(-10, -2), direction = 1)

# Save DET isoform enrichment plots and results
pdf(file.path(output_folder, "DET_gene_GO_enrichment_network.pdf"), width = 10, height = 10)
print(p5)
print(p6)
print(p7)
print(p8)
dev.off()

write.csv(ego3@result, file.path(output_folder, "DET_gene_GO_enrichment_network.csv"), row.names = FALSE)
