################### Di Wu ####################
################# July 2019 ###################

library(Seurat)
library(ggplot2)
library(Matrix)
#setwd("~/Documents/R_analysis/scRNA/EK-6744--04--17--2019/leukemia")

args=commandArgs(TRUE)
path<-args[1]
name<-args[2]
number<-args[3]

#matrix_dir = "/Users/wud3/Documents/R_analysis/scRNA/EK-6744--04--17--2019/leukemia/filtered_feature_bc_matrix/"
matrix_dir = path
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

pattern <- c("Hashtag1","Hashtag2","Hashtag3","Hashtag4","Hashtag5", "Hashtag6")
mat_tag <- mat[grep(paste(pattern, collapse="|"), rownames(mat), value = TRUE), ]
mat_RNA <- mat[1:27998, ]



# Confirm that the HTO have the correct names
#rownames(pbmc.htos)
rownames(mat_tag)
# Setup Seurat object
data.umis <- mat_RNA
data.hashtag <- mat_tag
data.hashtag <- CreateSeuratObject(counts = data.umis)

# Normalize RNA data with log normalization
data.hashtag <- NormalizeData(data.hashtag)
# Find and scale variable features
data.hashtag <- FindVariableFeatures(data.hashtag, selection.method = "mean.var.plot")
data.hashtag <- ScaleData(data.hashtag, features = VariableFeatures(data.hashtag))

# Add HTO data as a new assay independent from RNA
data.hashtag[["HTO"]] <- CreateAssayObject(counts = mat_tag)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
data.hashtag <- NormalizeData(data.hashtag, assay = "HTO", normalization.method = "CLR")
#Demultiplex cells based on HTO enrichment
data.hashtag <- HTODemux(data.hashtag, assay = "HTO", positive.quantile = 0.99)
# Global classification results
table(data.hashtag$HTO_classification.global)

Idents(data.hashtag) <- "HTO_maxID"
RidgePlot(data.hashtag, assay = "HTO", features = rownames(data.hashtag[["HTO"]])[1:2], ncol = 2)
FeatureScatter(data.hashtag, feature1 = "hto_Hashtag1", feature2 = "hto_Hashtag2")

#Compare number of UMIs for singlets, doublets and negative cells

Idents(data.hashtag) <- "HTO_classification.global"
VlnPlot(data.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
#VlnPlot(data.hashtag, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)
#Generate a two dimensional tSNE embedding for HTOs.Here we are grouping cells by singlets and doublets for simplicity.
# First, we will remove negative cells from the object
data.hashtag.subset <- subset(data.hashtag, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = data.hashtag.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
dara.hashtag.subset <- RunTSNE(data.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
#DimPlot(data.hashtag.subset)
HTOHeatmap(data.hashtag, assay = "HTO", ncells = 5000)
HTOHeatmap(data.hashtag, assay = "HTO")
write.csv(data.hashtag@meta.data, "hashtag_metadata_p0.99.csv")
hashtag.raw.data <- as.matrix(GetAssayData(data.hashtag@assays$HTO, slot = "counts"))
write.csv(t(hashtag.raw.data), "hashtag_raw_data.csv")

# Use hashtag_raw_data.csv to double check negative and then update metadata
a1 <- read.csv("hashtag_metadata_p0.99.csv")
head(a1)
data.hashtag@meta.data$hash.ID <- a1$hash.ID
data.hashtag@meta.data$HTO_classification.global <- a1$HTO_classification.global
data.hashtag@meta.data$HTO_classification <- a1$HTO_classification

pdf(paste("HTO","_heatmap_scRNA.pdf",sep=""),16,10, title = "test")
HTOHeatmap(data.hashtag, assay = "HTO", ncells = 5000)
dev.off()
