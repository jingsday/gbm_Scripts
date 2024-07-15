library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
# get data location
dirs <- list.files(path = '~/project_GBM/gbm_DATA/gbm_DATA_single_atlas', recursive = F, full.names = F)
dirs <-grep('barcodes.tsv.gz$',dirs,value=TRUE)
dirs
names_list <- c()  # Create an empty list to store the names

#names for reading documents
for (x in dirs) {
  name <- unique(sub('_barcodes.tsv.gz', '',x))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}


#names_list <-unique(names_list)
names_list <- names_list[names_list != "GSM5319529_SF8963"]
names_list
# Add "GSM5288674_Sample-11_" to the list
#names_list <- c(names_list, "GSM5288674_Sample-11_")

for (name in names_list) {
  cts <- ReadMtx(mtx = paste0(name,'_matrix.mtx.gz'),
                 features = paste0(name,'_features.tsv.gz'),
                 cells = paste0(name,'_barcodes.tsv.gz'),feature.column = 1)
  
  # create seurat objects
  assign(str_sub(name,start=12), CreateSeuratObject(counts = cts))
  
}
#str_sub(name,start=12), CreateSeuratObject(counts = cts))

ls()


merged_seurat <- merge(
  x = SF10099,
y=c(SF10099v2,  SF10108,SF10432,SF10433,SF10433v2,  
SF10441,SF10484,SF10565,SF10565v2, SF10592,SF10857,SF11082,SF11248,SF11344, SF11488,SF11587,SF11720,SF11720v2,  
SF11780,SF11815,SF11857,SF11873,SF11916,SF11981, SF1199, SF12008,
SF12115,SF12165,SF12165v2, SF12243,SF12243v2,  
SF12382,SF12407,SF12408, SF12460,SF12751,SF12754,SF1343, 
SF2501, SF2628, SF2777, SF2979, SF2990, SF3073, SF3076, 
SF3243, SF3391, SF3448, SF3996,SF4209, SF4209v2,SF4324, SF4449, SF4449v2,  
SF4810, SF4810v2,SF4849, SF5581, SF6098,SF6118, 
SF6118v2,SF6186, SF6621, SF6809, SF7025, SF7062, SF7307, SF7307v2,SF7388,
SF9358, SF9372, SF9494, SF9510, SF9715, SF9715v2,SF9798, SF9871, SF9962),
add.cell.ids=ls()[6:84]
,project='GBM')


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
saveRDS(seurat.integrated,'/project_GBM/gbm_OUTPUT/gbm_OUTPUT_paired_merged.rds')
