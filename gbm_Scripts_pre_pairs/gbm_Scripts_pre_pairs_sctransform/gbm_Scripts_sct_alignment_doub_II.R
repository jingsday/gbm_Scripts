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
obj.list <- SplitObject(merged_seurat, split.by = 'Sample')

#Part II: SCT intergration and save files for furthur analysis

obj.list <- lapply(X = obj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

seurat.integrated<- readRDS('~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_intergration_doubletII.rds')
saveRDS(seurat.integrated,file='~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_intergration_doubletII.rds')
########################################################################################

seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, verbose = FALSE,resolution = 0.5)
DimPlot(seurat.integrated, group.by = 'Sample')

head(seurat.integrated@meta.data)

tumor_info<-read.table("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_metadata/GSE174554_Tumor_normal_metadata_11916v2.txt", header = TRUE, sep = " ")
head(tumor_info)
head(seurat.integrated@meta.data)

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
seurat.integrated

saveRDS(seurat.integrated,file='~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_intergration_doublet_full.rds')

DefaultAssay(seurat.integrated) <- "SCT"
plots <-VlnPlot(seurat.integrated, features = c("PTPRZ1", 'VEGFA','SLC44A1' ), split.by = "Condition",group.by = 'Condition',
                pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 1)

DimPlot(seurat.integrated, reduction = "umap",group.by = 'Sample')
DimPlot(seurat.integrated, reduction = "umap",group.by = 'Condition')
DimPlot(seurat.integrated, reduction = "umap",group.by = 'Tumor_annotation')
DimPlot(seurat.integrated, reduction = "umap",group.by ='seurat_clusters')


###remove microglias based on biomakers expression
Idents(seurat.integrated) <- "Condition"

DefaultAssay(seurat.integrated) <- "SCT"
table(Idents(seurat.integrated))
seurat.integrated <- PrepSCTFindMarkers(seurat.integrated)

tmz.response <- FindMarkers(seurat.integrated, assay = "SCT", ident.1 = "Primary", ident.2 = "Recurrent", verbose = FALSE)
head(tmz.response, n = 15)

VlnPlot(seurat.integrated, features = c("PTPRZ1", "VEGFA", "P2RY12"), split.by = 'Condition',
        group.by = "Tumor_annotation", pt.size = 0, combine = FALSE)


#DEGs recommended to use RNA assay 
DefaultAssay(seurat.integrated) <- "integrated"


# Run the standard workflow for visualization and clustering
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:15)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
#sct_murine <- FindClusters(sct_murine, resolution = 0.5)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)

seurat.integrated.markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE)
seurat.integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.csv(seurat.integrated.markers, "seurat.integrated.markers_05.csv", row.names = FALSE)
seurat.integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(seurat.integrated, features = top10$gene) + NoLegend()

dev.off()
##downsize
# Downsample the data (if your dataset is large)
test <- subset(seurat.integrated, cells = sample(Cells(seurat.integrated), 3000))

DoHeatmap(test, features = top10$gene) + NoLegend()

table(seurat.integrated@meta.data$seurat_clusters)

DimPlot(seurat.integrated, reduction = "umap",group.by = 'Sample')
DimPlot(seurat.integrated, reduction = "umap",group.by = 'Condition')
DimPlot(seurat.integrated, reduction = "umap",group.by = 'Tumor_annotation')
DimPlot(seurat.integrated, reduction = "umap",group.by ='seurat_clusters')

cluster4.markers <- FindMarkers(seurat.integrated, ident.1 = 4)
head(cluster4.markers, n = 5)

table(seurat.integrated@meta.data$seurat_clusters)
  ============================
Idents(seurat.integrated) <- "seurat_clusters"
library(metap)
seurat.integrated.markers <- FindConservedMarkers(seurat.integrated, verbose= False)

seurat.integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

seurat.integrated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat.integrated, features = top10$gene) + NoLegend()
