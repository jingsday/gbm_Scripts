library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(ggplot2)
library(patchwork)


seurat.integrated <- readRDS('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_intergration.rds')
seurat.integrated

tumor_info<-read.table("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_metadata/GSE174554_Tumor_normal_metadata_11916v2.txt", header = TRUE, sep = " ")
head(tumor_info)
head(seurat.integrated@meta.data)

tumor_info$Barcode <- paste0(tumor_info$Barcode,'-','1')

merged_metadata <- seurat.integrated@meta.data %>%
  left_join(tumor_info, by = c("Sample", "Barcode"))
seurat.integrated@meta.data$Tumor_annotation <-merged_metadata$Tumor_Normal_annotation

head(seurat.integrated@meta.data)
seurat.integrated@meta.data$Condition <-'Recurrent'

seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF2777'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF2990'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF3076'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF3391'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF11916'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF11082'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF9358'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF9798'] <- 'Primary'


seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.1)

DimPlot(seurat.integrated, reduction = "umap",group.by = 'Sample')
DimPlot(seurat.integrated, reduction = "umap",split.by = 'Condition')

DimPlot(seurat.integrated, reduction = "umap",group.by = 'Condition')

DefaultAssay(seurat.integrated) <- "SCT"
plots <-VlnPlot(seurat.integrated, features = c("PTPRZ1", 'VEGFA','SLC44A1' ), split.by = "Condition",group.by = 'Condition',
                 pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 1)
