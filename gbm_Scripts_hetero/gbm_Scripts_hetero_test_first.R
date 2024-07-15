library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

setwd('')

#locations of files
mtx <- "/Users/lidiayung/project/resource/GSE174554_RAW/GSM5319548_SF2979_matrix.mtx.gz"
cells <- "/Users/lidiayung/project/resource/GSE174554_RAW/GSM5319548_SF2979_barcodes.tsv.gz"
features <- "/Users/lidiayung/project/resource/GSE174554_RAW/GSM5319548_SF2979_features.tsv.gz"
str(features)
sf2979 <- ReadMtx(mtx = mtx, cells = cells, features = features,feature.column = 1)


# Initialize the Seurat object with the raw (non-normalized data).
sf2979_p <- CreateSeuratObject(counts = sf2979, project = "sf2979", min.cells = 3, min.features = 200)
str(sf2979_p)


#read metadata and creating tumor_vector to reference
metadata <- read.csv("/Users/lidiayung/project/resource/GSE174554_RAW/GSE174554_Tumor_normal_metadata.txt", header = TRUE,sep=" ")
metadata_subset_sf2979 <- metadata[metadata$Sample == "SF2979", ]
metadata_subset <-metadata_subset_sf2979[metadata_subset_sf2979$Tumor_Normal_annotation == "Tumor", ]
head(metadata_subset)
str(metadata_subset)

tumor_vector<- paste0(metadata_subset$Barcode, "-1")
str(tumor_vector)

#dataset prepared to use
tumor_sf2979 <- subset(sf2979_p, cells = tumor_vector)
tumor_sf2979

tumor_sf2979[["percent.mt"]] <- PercentageFeatureSet(tumor_sf2979, pattern = "^MT-")
