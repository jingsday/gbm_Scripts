library(Seurat)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(Matrix)
library(stringr)
library(DoubletFinder)
library(tidyverse)

library(tidyverse)
library(DoubletFinder)

wkdir <- '~/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas/'

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
for (obj_name in c("SF11082", "SF11488","SF11082","SF11488","SF11916","SF12382","SF2777"
,"SF2979","SF2990","SF3073","SF3076","SF3243","SF3391","SF3448","SF9358","SF9494","SF9798","SF9962")) {

  obj <- get(obj_name)  # Retrieve the Seurat object using its name

  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

  obj <- NormalizeData(object = obj)
  obj <- FindVariableFeatures(object = obj)
  obj <- ScaleData(object = obj)
  obj <- RunPCA(object = obj)
  ElbowPlot(obj)
  obj <- FindNeighbors(object = obj, dims = 1:20)
  obj <- FindClusters(object = obj)
  obj <- RunUMAP(object = obj, dims = 1:20)

  ## pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)

  pK <- as.numeric(as.character(pK[[1]]))
  print(pK)

  ## Homotypic Doublet Proportion Estimate
  annotations <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.015 * nrow(obj@meta.data))  ## Assuming ~1.5% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  ## Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

  ## Visualize doublets
  p1 <- DimPlot(obj, reduction = 'umap', group.by = colnames(obj@meta.data)[ncol(obj@meta.data)])
  
  # Save UMAP plot for doublet detection
  ggsave(plot = p1, filename = paste0('~/Phd_project/project_GBM/gbm_Scripts/gbm_Scripts_pre_pairs/gbm_Scripts_pre_pairs_OUTPUT/', obj_name, '_doublet.png'))

  ## Output metadata and table of doublets
  meta_col <- colnames(obj@meta.data)[ncol(obj@meta.data)]
  paste(str(obj), table(obj@meta.data[meta_col]))

  # Save the modified Seurat object back into the original name
  assign(obj_name, obj)
}

###Retain singlet
for (obj in c(SF11082,SF11488,SF11916,SF12382,SF2777,SF2979,SF2990,SF3073,SF3076,SF3243,SF3391,SF3448,SF9358,SF9494,SF9798,SF9962)){

  Idents(obj) <- colnames(obj@meta.data)[ncol(obj@meta.data)]
  print(table(Idents(obj)))
  obj <- subset(x = obj, idents = "Singlet")
}


for (obj in c(SF11082,SF11488,SF11916,SF12382,SF2777,SF2979,SF2990,SF3073,SF3076,SF3243,SF3391,SF3448,SF9358,SF9494,SF9798,SF9962)){
  obj_name
  obj$percent.mt <- PercentageFeatureSet(obj,pattern = "^MT-")

  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

  obj <- NormalizeData(object = obj)
  obj <- FindVariableFeatures(object = obj)
  obj <- ScaleData(object = obj)
  obj <- RunPCA(object = obj)
  ElbowPlot(obj)
  obj <- FindNeighbors(object = obj, dims = 1:20)
  obj <- FindClusters(object = obj)
  obj <- RunUMAP(object = obj, dims = 1:20)
## pK Identification (no ground-truth) ------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 

  pK <- as.numeric(as.character(pK[[1]]))
  print(pK)
## Homotypic Doublet Proportion Estimate ------------------------------------------

  annotations <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.015*nrow(obj@meta.data))## Assuming ~% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  nExp_poi.adj
# run doubletFinder 
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)
                               
  head(obj@meta.data)                        
# visualize doublets
  p1 <- DimPlot(obj, reduction = 'umap', group.by = colnames(obj@meta.data)[ncol(obj@meta.data)])
  
  ggsave(plot=p1,filename = '~/Phd_project/project_GBM/gbm_Scripts/doublet.png')

  meta_col <-colnames(obj@meta.data)[ncol(obj@meta.data)]

  paste(str(obj),table(obj@meta.data[meta_col]))
}





#Doublet 
SF11082.filtered <- NormalizeData(object = SF11082)
SF11082.filtered <- FindVariableFeatures(object = SF11082.filtered)
SF11082.filtered <- ScaleData(object = SF11082.filtered)
SF11082.filtered <- RunPCA(object = SF11082.filtered)
ElbowPlot(SF11082.filtered)
SF11082.filtered <- FindNeighbors(object = SF11082.filtered, dims = 1:20)
SF11082.filtered <- FindClusters(object = SF11082.filtered)
SF11082.filtered <- RunUMAP(object = SF11082.filtered, dims = 1:20)



## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_SF11082 <- paramSweep(SF11082.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_SF11082 <- summarizeSweep(sweep.res.list_SF11082, GT = FALSE)
bcmvn_SF11082 <- find.pK(sweep.stats_SF11082)

p1 <- ggplot(bcmvn_SF11082, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

outdir <- '~/Phd_project/project_GBM/gbm_Scripts/gbm_Scripts_pre_pairs/gbm_Scripts_pre_pairs_sctransform/00_output_doublet/'
ggsave(filename=file.path(outdir,'doublet_SF11082.png'),plot=p1,width=15)

pK <- bcmvn_SF11082 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- SF11082.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.015*nrow(SF11082.filtered@meta.data))  ## Assuming 1.5% doublet formation rate =
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
SF11082.seurat.filtered <- doubletFinder_v3(SF11082.seurat.filtered, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(SF11082.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")


# number of singlets and doublets
table(SF11082.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)




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
  add.cell.ids = ls()[15:30],project = 'GBM')


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
ggsave('features.jpg')


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

saveRDS(seurat.integrated,file='~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_intergration_doublet.rds')


seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, verbose = FALSE)
p2 <- DimPlot(seurat.integrated, group.by = 'Sample')
ggsave(plot=p2,filename = '~/Phd_project/project_GBM/gbm_Scripts/integration_sct_origin_doublet.png')
head(seurat.integrated@meta.data)
