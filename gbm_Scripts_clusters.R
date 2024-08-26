library(Seurat)
library(dplyr)
data_dir <- '/Users/lidiayung/PhD_project/project_GBM/gbm_DATA/sct_cl00.rds'

int_mouse <- readRDS(data_dir)
head(int_mouse@meta.data)

head(Idents(int_mouse))
Idents(int_mouse) <- int_mouse@meta.data$Sample

GSM5288668 <- subset(x = int_mouse, idents = "GSM5288668")
GSM5288668_down <- subset(GSM5288668_down, subset = nFeature_RNA > 200 & nFeature_RNA < 8000)
GSM5288668_down <- subset(GSM5288668_down, subset = nCount_RNA > 500)
sum(is.na(GSM5288668_down))



GSM5288668_down <- RunPCA(GSM5288668, verbose = FALSE)
GSM5288668_down <- RunUMAP(GSM5288668_down, dims = 1:30, verbose = FALSE)

GSM5288668_down <- FindNeighbors(GSM5288668_down, dims = 1:30, verbose = FALSE)
GSM5288668_down <- FindClusters(GSM5288668_down, ,resolution=0.1)
DimPlot(GSM5288668_down, label = TRUE)

table(GSM5288668_down$seurat_clusters)



GSM5288668_down.markers <- FindAllMarkers(GSM5288668_down, only.pos = FALSE)
GSM5288668_down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(abs(avg_log2FC > 1))

GSM5288668_down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(GSM5288668_down, features = top10$gene) + NoLegend()
hist(GSM5288668_down.markers$avg_log2FC)
VlnPlot(GSM5288668_down, features = c("Krt20", "Upk1b"),layer = 'SCT')



