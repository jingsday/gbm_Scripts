library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)


setwd('/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas')

names_list <- c('GSM5319518_SF2777','GSM5319548_SF2979','GSM5319519_SF2990',
                'GSM5319549_SF3073','GSM5319520_SF3076','GSM5319550_SF3243',
                'GSM5319521_SF3391','GSM5319551_SF3448','GSM5319511_SF11916',
                'GSM5319543_SF12382','GSM5319506_SF11082','GSM5319562_SF11488',
                'GSM5319530_SF9358','GSM5319568_SF9962','GSM5319559_SF9798','GSM5319532_SF9494')

for (name in names_list) {
  cts <- ReadMtx(mtx = paste0(name,'_matrix.mtx.gz'),
                 features = paste0(name,'_features.tsv.gz'),
                 cells = paste0(name,'_barcodes.tsv.gz'),feature.column = 1)
  
  # create seurat objects
  assign(str_sub(name,start=12), CreateSeuratObject(counts = cts))}
  
  
ls()

merged_seurat <- merge(
  x = SF11082,
  y = c(
        SF11488,SF11916,SF12382,SF2777,           
        SF2979,SF2990,SF3073,SF3076,SF3243,
         SF3391,SF3448,SF9358,SF9494,SF9798,
        SF9962),
  add.cell.ids = ls()[6:21],
  project = 'GBM'
)

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)
head(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')


obj.list <- SplitObject(merged_seurat, split.by = 'Sample')


obj.list <- lapply(X = obj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

#Data ready to be retrieved here 
saveRDS(seurat.integrated,"~/project_GBM/gbm_OUTPUT/gbm_project_integrated.rds") 