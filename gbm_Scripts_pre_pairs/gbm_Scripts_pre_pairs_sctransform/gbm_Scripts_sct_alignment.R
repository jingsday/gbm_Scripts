library(Seurat)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(Matrix)
library(stringr)

wkdir <- '/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas/'

names_list <- c('GSM5319518_SF2777','GSM5319548_SF2979','GSM5319519_SF2990',
                'GSM5319549_SF3073','GSM5319520_SF3076','GSM5319550_SF3243',
                'GSM5319521_SF3391','GSM5319551_SF3448','GSM5319511_SF11916',
                'GSM5319543_SF12382','GSM5319506_SF11082','GSM5319562_SF11488',
                'GSM5319530_SF9358','GSM5319568_SF9962','GSM5319559_SF9798','GSM5319532_SF9494')

for (name in names_list) {
  cts <- ReadMtx(mtx = paste0(wkdir,name,'_matrix.mtx.gz'),
                 features = paste0(wkdir,name,'_features.tsv.gz'),
                 cells = paste0(wkdir,name,'_barcodes.tsv.gz'),feature.column = 1)
  
  # create seurat objects
  assign(str_sub(name,start=12), CreateSeuratObject(counts = cts,min.cells = 3, min.features = 200))
  
}

ls()

# Create list of gene names from all Seurat objects
total.genes <- list(rownames(SF11082),
                    rownames(SF11488),
                    rownames(SF11916),
                    rownames(SF12382),
                    rownames(SF2777),
                    rownames(SF2979),
                    rownames(SF2990),
                    rownames(SF3073),
                    rownames(SF3076),
                    rownames(SF3243),
                    rownames(SF3391),
                    rownames(SF3448),
                    rownames(SF9358),
                    rownames(SF9494),
                    rownames(SF9798),
                    rownames(SF9962))

# Get common gene names 
common.genes <- Reduce(f = intersect, x = total.genes)

SF11082 <- SF11082[common.genes,]
SF11488 <- SF11488[common.genes,]
SF11916 <- SF11916[common.genes,]
SF12382 <- SF12382[common.genes,]
SF2777 <- SF2777[common.genes,]
SF2979 <- SF2979[common.genes,]
SF2990 <- SF2990[common.genes,]
SF3073 <- SF3073[common.genes,]
SF3076 <- SF3076[common.genes,]
SF3243 <- SF3243[common.genes,]
SF3391 <- SF3391[common.genes,]
SF3448 <- SF3448[common.genes,]
SF9358 <- SF9358[common.genes,]
SF9494 <- SF9494[common.genes,]
SF9798 <- SF9798[common.genes,]
SF9962 <- SF9962[common.genes,]

ls()

merged_seurat <- merge(
  x = SF11082,y = c(SF11488,SF11916,SF12382,SF2777,
                   SF2979,SF2990,SF3073,SF3076,SF3243,
                   SF3391,SF3448,SF9358,SF9494,SF9798,SF9962),
  add.cell.ids = ls()[5:20],project = 'GBM')


head(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)
head(merged_seurat@meta.data)


# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^Mt-')

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#ggsave('features.jpg')


merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                                   nFeature_RNA < 5000 &
                                   mitoPercent < 2.5)
merged_seurat
merged_seurat_filtered


obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Sample')

#Part II: SCT intergration and save files for furthur analysis

obj.list <- lapply(X = obj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(seurat.integrated,file='~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_intergration.rds')


seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, verbose = FALSE)
DimPlot(seurat.integrated, label = TRUE)
head(seurat.integrated@meta.data)
