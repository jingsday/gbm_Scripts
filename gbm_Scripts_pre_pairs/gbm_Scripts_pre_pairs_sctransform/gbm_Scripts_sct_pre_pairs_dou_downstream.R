library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(patchwork)

dou_seurat.integrated <- readRDS('~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_intergration_doublet.rds')
dou_seurat.integrated

tumor_info<-read.table("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_metadata/GSE174554_Tumor_normal_metadata_11916v2.txt", header = TRUE, sep = " ")
head(tumor_info)
head(dou_seurat.integrated@meta.data)

tumor_info$Barcode <- paste0(tumor_info$Barcode,'-','1')

merged_metadata <- dou_seurat.integrated@meta.data %>%
  left_join(tumor_info, by = c("Sample", "Barcode"))
dou_seurat.integrated@meta.data$Tumor_annotation <-merged_metadata$Tumor_Normal_annotation


dou_seurat.integrated@meta.data$Condition <-'Recurrent'

dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF2777'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF2990'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF3076'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF3391'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF11916'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF11082'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF9358'] <- 'Primary'
dou_seurat.integrated@meta.data$Condition[dou_seurat.integrated@meta.data$Sample == 'SF9798'] <- 'Primary'
dou_seurat.integrated

#SCT assay
DefaultAssay(dou_seurat.integrated) <- "SCT"

dou_seurat.integrated <- ScaleData(object = dou_seurat.integrated)
dou_seurat.integrated <- RunPCA(object = dou_seurat.integrated)
dou_seurat.integrated <- RunUMAP(object = dou_seurat.integrated, dims = 1:50)

dou_seurat.integrated <- FindNeighbors(dou_seurat.integrated, dims = 1:30, verbose = FALSE)
dou_seurat.integrated <- FindClusters(dou_seurat.integrated, verbose = FALSE,resolution = 0.1)

plots <-VlnPlot(dou_seurat.integrated, features = c("PTPRZ1", 'VEGFA','SLC44A1' ), split.by = "Condition",group.by = 'Condition',
                pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 1)

DimPlot(dou_seurat.integrated, reduction = "umap",group.by = 'Sample')
DimPlot(dou_seurat.integrated, reduction = "umap",group.by = 'Condition')
DimPlot(dou_seurat.integrated, reduction = "umap",group.by = 'Tumor_annotation')
DimPlot(dou_seurat.integrated, reduction = "umap",group.by ='seurat_clusters')


DefaultAssay(dou_seurat.integrated) <- "integrated"

dou_seurat.integrated <- ScaleData(object = dou_seurat.integrated)
dou_seurat.integrated <- RunPCA(object = dou_seurat.integrated)
dou_seurat.integrated <- RunUMAP(object = dou_seurat.integrated, dims = 1:50)

dou_seurat.integrated <- FindNeighbors(dou_seurat.integrated, dims = 1:30, verbose = FALSE)
dou_seurat.integrated <- FindClusters(dou_seurat.integrated, verbose = FALSE,resolution = 0.1)
DimPlot(dou_seurat.integrated, group.by = 'Sample')


#DEGs recommended to use RNA assay 
DefaultAssay(dou_seurat.integrated) <- "integrated"


dou_seurat.integrated <- ScaleData(dou_seurat.integrated, verbose = FALSE)

dou_seurat.integrated <- RunPCA(dou_seurat.integrated, verbose = FALSE)
dou_seurat.integrated <- RunUMAP(dou_seurat.integrated, reduction = "pca", dims = 1:30)

dou_seurat.integrated <- FindNeighbors(dou_seurat.integrated, dims = 1:30, verbose = FALSE)

dou_seurat.integrated <- FindClusters(dou_seurat.integrated, resolution = 0.1)

dou_seurat.integrated.markers <- FindAllMarkers(dou_seurat.integrated, only.pos = TRUE)
dou_seurat.integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

dou_seurat.integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
library(dplyr)
# Ensure correct assay and data slot
expr_matrix <- GetAssayData(dou_seurat.integrated, assay = "RNA", slot = "data")

DoHeatmap(dou_seurat.integrated, features = top10$gene) 
DoHeatmap(dou_seurat.integrated, features = top10$gene, size = 3) + NoLegend()

top10[top10$cluster==4,]
FeaturePlot(dou_seurat.integrated, features = c("CD44", "SNTG1", "SLC11A1", "CLDN11", "DIAPH3", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
VlnPlot(dou_seurat.integrated, features = top10$gene[2])
install.packages(c('pheatmap','ComplexHeatmap'))

library(ComplexHeatmap)
heatmap_data <- dou_seurat.integrated.markers %>%
  select(gene, cluster, avg_log2FC) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene")

# Create heatmap using pheatmap
pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE)
Heatmap(dou_seurat.integrated.markers[c('avg_log2FC','gene'),], celltype=dou_seurat.integrated.markers$cluster)
length(dou_seurat.integrated.markers$cluster)
# Create the heatmap
pheatmap(subset_matrix, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = FALSE)

png('test.png')
DoHeatmap(dou_seurat.integrated, features = top10$gene) + NoLegend()
dev.off()
dou_seurat.integrated
#Singlet retrieving


singlets_dou_seurat.integrated <- RunPCA(singlets_dou_seurat.integrated, verbose = FALSE)
singlets_dou_seurat.integrated <- RunUMAP(singlets_dou_seurat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
singlets_dou_seurat.integrated <- FindNeighbors(singlets_dou_seurat.integrated, reduction = "pca", dims = 1:30)
singlets_dou_seurat.integrated <- FindClusters(singlets_dou_seurat.integrated, resolution = 0.1)


###remove microglias based on biomakers expression
Idents(singlets_dou_seurat.integrated) <- "Condition"

DefaultAssay(singlets_dou_seurat.integrated) <- "SCT"
table(Idents(singlets_dou_seurat.integrated))
singlets_dou_seurat.integrated <- PrepSCTFindMarkers(singlets_dou_seurat.integrated)

library(future)
options(future.globals.maxSize = 10 * 1024^3) 

tmz.response <- FindMarkers(singlets_dou_seurat.integrated, assay = "SCT", ident.1 = "Primary", ident.2 = "Recurrent", verbose = FALSE, recorrect_umi = FALSE)

head(tmz.response, n = 15)

VlnPlot(dou_seurat.integrated, features = c("PTPRZ1", "VEGFA", "SLC44A1"), split.by = 'Condition',
                 group.by = "Tumor_annotation", pt.size = 0, combine = FALSE)

