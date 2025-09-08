#!/usr/bin/env Rscript
# Script: visualize novel isoforms in TG and their existence in other tissues.r
# Description: Visualizes the expression levels of novel isoforms of TG in TG and mouse encode longread RNA-seq tissues using t-SNE and ComplexHeatmap.

## 1. Load Required Libraries 
library(Seurat)
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggrepel)

## 2. User-defined Paths (edit as needed) 
data_file <- "data/longread_novel_isoform_counts.tsv"# file containing the expression data of novel isoforms of TG and mouse encode longread RNA-seq tissues
sample_info_file <- "data/mouse_longread_sample_info.tsv"# file containing sample information of mouse encode longread RNA-seq tissues
output_dir <- "results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## 3. Load Data 
expr_all <- read.table(data_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(expr_all) <- expr_all[, 1]

# Extract counts matrix
sample_cols <- 2:(ncol(expr_all) - 3)
counts <- as.matrix(expr_all[, sample_cols])
storage.mode(counts) <- "integer"

## 4. Create Seurat Object 
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

## 5. t-SNE Visualization 
tsne_plot <- DimPlot(
  seurat_obj, reduction = "tsne", group.by = "orig.ident",
  label = FALSE, pt.size = 1
) +
  ggtitle("t-SNE of Samples by Tissue Group") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "longread_isoform_TSNE.pdf"), tsne_plot, width = 5, height = 4)

## 6. Heatmap Visualization 
scaled_mat <- as.matrix(seurat_obj@assays$RNA$scale.data)

mean_expr_mat <- do.call(cbind, lapply(unique(seurat_obj$orig.ident), function(grp) {
  samples <- colnames(seurat_obj)[seurat_obj$orig.ident == grp]
  rowMeans(scaled_mat[, samples, drop = FALSE])
}))

# K-means clustering
set.seed(123)
num_clusters <- 10
kmeans_result <- kmeans(mean_expr_mat, centers = num_clusters)
row_split <- kmeans_result$cluster

ht <- Heatmap(
  mean_expr_mat,
  name = "Z-score",
  col = colorRamp2(seq(-2, 2, length.out = 50), colorRampPalette(c("#4DBBD5FF", "white", "#ED0000FF"))(50)),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  row_split = row_split
)

pdf(file.path(output_dir, "longread_isoform_Heatmap.pdf"), width = 5, height = 10)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

## 7. Save Isoform Cluster Assignments 
isoform_clusters_df <- data.frame(
  cluster = row_split,
  isoform = names(row_split)
)
write.table(
  isoform_clusters_df,
  file = file.path(output_dir, "novel_isoform_clusters.tsv"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)

message("Analysis complete. Results saved to: ", output_dir)
