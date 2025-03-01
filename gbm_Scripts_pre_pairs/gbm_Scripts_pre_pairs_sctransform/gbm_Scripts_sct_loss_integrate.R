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
  add.cell.ids = ls()[4:19],
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
saveRDS(seurat.integrated,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_project_integrated.rds") 
head(seurat.integrated@meta.data)

seu <- subset(seurat.integrated,subset = Sample == 'SF11082')  

#Primary and recurrent

metadata <-read.delim2('/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_metadata/GSE174554_Tumor_normal_metadata_11916v2.txt',
                      sep = " ")
head(metadata)
'SF11082_AAACCCAAGCATCAAA-1  '
metadata$Barcode <- paste0(metadata$Sample, "_", metadata$Barcode, "-1")
rownames(metadata) <- metadata$Barcode
head(rownames(metadata))

df <- metadata[metadata$Barcode %in% rownames(seurat.integrated@meta.data), ]
rownames(df) <- df$Barcode
seurat.integrated@meta.data$Tumor_Normal_annotation <- df[rownames(seurat.integrated@meta.data), "Tumor_Normal_annotation"]

seurat.integrated@meta.data$Status <- 'Recurrent'

# List of samples to change the 'Status' to 'Primary'
primary_samples <- c('SF2770', 'SF2990', 'SF3076', 'SF3391', 'SF11916', 'SF11082', 'SF9358', 'SF9798')

# Change 'Status' to 'Primary' for matching samples
for (sample in primary_samples) {
  seurat.integrated@meta.data$Status[grepl(sample, seurat.integrated@meta.data$Sample, fixed = TRUE)] <- 'Primary'
}

head(seurat.integrated@meta.data)

# Retrieve df for Support vector machine 

tumor_seu <- subset(seurat.integrated, subset = Tumor_Normal_annotation=='Tumor')
tumor_seu_prm <- subset(tumor_seu, subset = Status=='Primary')
tumor_seu_rct <- subset(tumor_seu, subset = Status=='Recurrent')

prm_corrected_UMI <- tumor_seu_prm[["SCT"]]$data
prm_corrected_UMI <-t(prm_corrected_UMI)

writeMM(prm_corrected_UMI, "/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_tumor_prm_corrected_UMI.mtx")
write.table(rownames(prm_corrected_UMI), 
            file = "/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_tumor_prm_corrected_UMI_cells.txt", col.names = FALSE, row.names = FALSE)
write.table(colnames(prm_corrected_UMI), 
            file = "/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_tumor_prm_corrected_UMI_genes.txt", col.names = FALSE, row.names = FALSE)


rct_corrected_UMI <- tumor_seu_rct[["SCT"]]$data
rct_corrected_UMI <-t(rct_corrected_UMI)

writeMM(rct_corrected_UMI, "/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_tumor_rct_corrected_UMI.mtx")
write.table(rownames(rct_corrected_UMI), 
            file = "/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_tumor_rct_corrected_UMI_cells.txt", col.names = FALSE, row.names = FALSE)
write.table(colnames(rct_corrected_UMI), 
            file = "/home/jing/Phd_project/project_GBM/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_tumor_rct_corrected_UMI_genes.txt", col.names = FALSE, row.names = FALSE)

DefaultAssay(seurat.integrated) <- "integrated"
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)


p1 <- DimPlot(seurat.integrated, reduction = "umap")
p2 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Sample", 
              repel = TRUE)

p1 + p2
