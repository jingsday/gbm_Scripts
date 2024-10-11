library(Seurat)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(Matrix)
library(stringr)
library(DoubletFinder)
library(tidyverse)
library(patchwork)

wkdir <- '~/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas/'

merged_seurat<- readRDS('~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_merged_seurat_doublet.rds')
merged_seurat

merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent.mt")

seurat.integrated <- SCTransform(merged_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(seurat.integrated,file='~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_intergration_v2.rds')
seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30, verbose = FALSE)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, verbose = FALSE,resolution = 0.05)
DimPlot(seurat.integrated, label = TRUE)  


tumor_info<-read.table("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_metadata/GSE174554_Tumor_normal_metadata_11916v2.txt", header = TRUE, sep = " ")
head(tumor_info)
head(seurat.integrated@meta.data)
#join tables
tumor_info$Barcode <- paste0(tumor_info$Barcode,'-','1')

merged_metadata <- seurat.integrated@meta.data %>%
  left_join(tumor_info, by = c("Sample", "Barcode"))
seurat.integrated@meta.data$Tumor_annotation <-merged_metadata$Tumor_Normal_annotation


seurat.integrated@meta.data$Condition <-'Recurrent'

seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF2777'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF2990'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF3076'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF3391'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF11916'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF11082'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF9358'] <- 'Primary'
seurat.integrated@meta.data$Condition[seurat.integrated@meta.data$Sample == 'SF9798'] <- 'Primary'


DimPlot(seurat.integrated, reduction = "umap",group.by = 'Sample')
DimPlot(seurat.integrated, reduction = "umap",group.by = 'Condition')
DimPlot(seurat.integrated, reduction = "umap",group.by = 'Tumor_annotation')
DimPlot(seurat.integrated, reduction = "umap",group.by ='seurat_clusters')


# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(seurat.integrated, features = c("EGFR",'MALAT1'), pt.size = 0.2,
            ncol = 3)
Idents(seurat.integrated) <- "Tumor_annotation"

DefaultAssay(seurat.integrated) <- "SCT"
table(Idents(seurat.integrated))
seurat.integrated <- PrepSCTFindMarkers(seurat.integrated)

tmz.response <- FindMarkers(seurat.integrated, assay = "SCT", ident.1 = "Tumor", ident.2 = "Normal", verbose = FALSE)
head(tmz.response, n = 15)

VlnPlot(seurat.integrated, features = c("PTPRZ1", "VEGFA", "P2RY12"), split.by = 'Tumor_annotation',
        group.by = "Tumor_annotation", pt.size = 0, combine = FALSE)
