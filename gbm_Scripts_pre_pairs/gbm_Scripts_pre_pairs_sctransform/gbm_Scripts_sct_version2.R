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

seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30, verbose = FALSE)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, verbose = FALSE)
DimPlot(seurat.integrated, label = TRUE)  