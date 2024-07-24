# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(sctransform)


setwd('/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas')

# get data location
dirs <- list.files(recursive = F, full.names = F)
dirs <-grep('barcodes.tsv.gz$',dirs,value=TRUE)


names_list <- c()  # Create an empty list to store the names
#names for reading documents
for (x in dirs) {
  name <- unique(sub('_barcodes.tsv.gz', '',x))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}


names_list <-c('GSM5319518_SF2777','GSM5319548_SF2979')

for (name in names_list) {
  cts <- ReadMtx(mtx = paste0(name,'_matrix.mtx.gz'),
                 features = paste0(name,'_features.tsv.gz'),
                 cells = paste0(name,'_barcodes.tsv.gz'),feature.column = 1)
  
  # create seurat objects
  assign(str_sub(name,start=12), CreateSeuratObject(counts = cts,min.cells = 3, min.features = 200))
  
}


#read metadata and creating tumor_vector to reference
metadata <- read.csv("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_metadata/GSE174554_Tumor_normal_metadata.txt", header = TRUE,sep=" ")
metadata_subset_sf2979 <- metadata[metadata$Sample == "SF2979", ]
metadata_subset <-metadata_subset_sf2979[metadata_subset_sf2979$Tumor_Normal_annotation == "Tumor", ]
head(metadata_subset)

str(metadata_subset)

tumor_vector<- paste0(metadata_subset$Barcode, "-1")
str(tumor_vector)#dataset prepared to use
tumor_SF2979 <- subset(SF2979, cells = tumor_vector)
tumor_SF2979
#tumor_SF2979[["percent.mt"]] <- PercentageFeatureSet(tumor_sf2979, pattern = "^MT-")


#read metadata and creating tumor_vector to reference

metadata_subset_sf2777 <- metadata[metadata$Sample == "SF2777", ]
metadata_subset <-metadata_subset_sf2777[metadata_subset_sf2777$Tumor_Normal_annotation == "Tumor", ]
head(metadata_subset)
str(metadata_subset)
tumor_sf2777 <- subset(SF2777, cells = metadata_subset$Barcode)
tumor_sf2777

ls()
merged_seurat <- merge(x = tumor_sf2777,y = tumor_SF2979,
                       add.cell.ids = ls()[12:13],project = 'GBM')

head(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)
table(merged_seurat@meta.data['Sample'])

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Tumor','Sample', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^Mt-')

#Murine cells were filtered to retain higher quality cells (>200 &<8000 uniquely identified genes
#<25% of reads mapped to mitochondrial genes),
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                                   mitoPercent < 5)

VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA",'mitoPercent'), ncol = 2)

obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Sample')

common_genes <- Reduce(intersect, lapply(obj.list, rownames))
length(common_genes)

obj.list_aligned <- lapply(obj.list, function(x) {
  x <- subset(x, features = common_genes)
  return(x)
})

merged_seurat_aligned <- merge(x = obj.list_aligned$SF2777,y = obj.list_aligned$SF2979,project = 'GBM')

obj <- NormalizeData(merged_seurat_aligned)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Sample"))
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
library(reticulate)
library(SeuratWrappers)
obj@assays
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/home/jing/miniforge3/envs/scvi", verbose = FALSE
)


merged_seurat_aligned <- lapply(X = merged_seurat_aligned, FUN = SCTransform)


options(future.globals.maxSize = 3e+09)
merged_seurat_aligned <- SCTransform(merged_seurat_aligned)
merged_seurat_aligned <- RunPCA(merged_seurat_aligned, npcs = 30, verbose = F)

features <- SelectIntegrationFeatures(object.list = merged_seurat_aligned, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = merged_seurat_aligned, anchor.features = features)



merged_seurat_aligned <- IntegrateLayers(
  object = merged_seurat_aligned,
  method = 
  normalization.method = "SCT",
  verbose = F
)
merged_seurat_aligned <- FindNeighbors(merged_seurat_aligned, dims = 1:30, reduction = "integrated.dr")
merged_seurat_aligned <- FindClusters(merged_seurat_aligned, resolution = 0.5,cluster.name = "sct_clusters")

merged_seurat_aligned <- RunUMAP(merged_seurat_aligned, dims = 1:30, reduction = "integrated.dr")

DimPlot(merged_seurat_aligned,reduction = 'integrated.dr')
DimPlot(merged_seurat_aligned,group.by=c('Sample','sct_clusters'))

obj@assays[['SCT']]
