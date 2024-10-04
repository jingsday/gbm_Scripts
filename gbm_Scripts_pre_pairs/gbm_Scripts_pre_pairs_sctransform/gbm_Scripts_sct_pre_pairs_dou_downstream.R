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




#Singlet retrieving


doublet_columns <- grep("DF", colnames(dou_seurat.integrated@meta.data), value = TRUE)
dou_seurat.integrated@meta.data$singlet_anno <- apply(dou_seurat.integrated@meta.data[, doublet_columns], 1, function(row) {
  # Return the first non-NA value in the row
  non_na_value <- row[!is.na(row)][1]
  if (is.na(non_na_value)) {
    return(NA)  # If all values are NA, keep it as NA
  } else {
    return(non_na_value)  # Otherwise, return the first non-NA value (Singlet or Doublet)
  }
})
table(dou_seurat.integrated@meta.data$singlet_anno)
#Subset
Idents(dou_seurat.integrated) <- dou_seurat.integrated@meta.data$singlet_anno

singlets_seurat.integrated <- subset(dou_seurat.integrated,idents = 'Singlet')
singlets_seurat.integrated
#
head(singlets_seurat.integrated@meta.data)

Idents(singlets_seurat.integrated) <- dou_seurat.integrated@meta.data$orig.ident


singlets_seurat.integrated <- RunPCA(singlets_seurat.integrated, verbose = FALSE)
singlets_seurat.integrated <- RunUMAP(singlets_seurat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
singlets_seurat.integrated <- FindNeighbors(singlets_seurat.integrated, reduction = "pca", dims = 1:30)
singlets_seurat.integrated <- FindClusters(singlets_seurat.integrated, resolution = 0.1)

DefaultAssay(singlets_seurat.integrated) <- "SCT"
plots <-VlnPlot(singlets_seurat.integrated, features = c("PTPRZ1", 'VEGFA','SLC44A1' ), split.by = "Condition",group.by = 'Condition',
                pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 1)

DimPlot(singlets_seurat.integrated, reduction = "umap",group.by = 'Sample')
DimPlot(singlets_seurat.integrated, reduction = "umap",group.by = 'Condition')
DimPlot(singlets_seurat.integrated, reduction = "umap",group.by = 'Tumor_annotation')
DimPlot(singlets_seurat.integrated, reduction = "umap",group.by ='seurat_clusters')

###remove microglias based on biomakers expression
Idents(singlets_seurat.integrated) <- "Condition"

DefaultAssay(singlets_seurat.integrated) <- "SCT"
table(Idents(singlets_seurat.integrated))
singlets_seurat.integrated <- PrepSCTFindMarkers(singlets_seurat.integrated)

library(future)
options(future.globals.maxSize = 10 * 1024^3) 

tmz.response <- FindMarkers(singlets_seurat.integrated, assay = "SCT", ident.1 = "Primary", ident.2 = "Recurrent", verbose = FALSE, recorrect_umi = FALSE)

head(tmz.response, n = 15)

VlnPlot(dou_seurat.integrated, features = c("PTPRZ1", "VEGFA", "SLC44A1"), split.by = 'Condition',
                 group.by = "Tumor_annotation", pt.size = 0, combine = FALSE)

