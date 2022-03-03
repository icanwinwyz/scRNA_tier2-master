################# Di Wu ####################
############## Feb, 2020 ##################

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(gplots)
packageVersion("Seurat")


args=commandArgs(TRUE)
path<-args[1]
name<-args[2]
number<-args[3]
species <- arg[4]

test.data <- Read10X(data.dir=path)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 0 cells . Keep all cells with at
# least 300 detected genes
test <- CreateSeuratObject(counts = test.data, min.cells = 0, project = name)
test

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^Mt-")
if (species == "human") {
	test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^MT-")
	} else if (species == "mouse") {
		test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^mt-")
	} else if (species == "rat") {
		test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^Mt-")
	}



# Show QC metrics for the first 5 cells
head(test@meta.data, 5)

# Visualize QC metrics as a violin plot
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1, plot2))

#We filter cells that have unique feature counts less than 300
#We filter cells that have >15% mitochondrial counts
test_1 <- subset(test, subset = nFeature_RNA >= 300 & percent.mt <= 15)

plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
test_4 <- subset(test, subset = nFeature_RNA < 300 & percent.mt > 15)
e<-dim(test_4@assays$RNA)[2]
test_2 <- subset(test, subset = percent.mt <= 15)
a<-dim(test@assays$RNA)[2]-dim(test_2@assays$RNA)[2]-e
test_3 <- subset(test, subset = nFeature_RNA >= 300)
b<-dim(test@assays$RNA)[2]-dim(test_3@assays$RNA)[2]-e
c<-dim(test_1@assays$RNA)[2]
d<-dim(test@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells failed mito% < 15%",sep=" ")
text3<-paste(b,"cells failed total # expressed genes > 300.",sep=" ")
text5<-paste(e,"cells failed total # expressed genes > 300 and mito% < 15%.",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");

textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()

temp<-test_1@assays$RNA
dim(temp)
#temp<-temp[,which(temp[grep(paste("^","PECAM1","$",sep=""),rownames(temp),ignore.case=T),]==0)]
#dim(temp)
#temp<-temp[,which(temp[grep(paste("^","PTPRC","$",sep=""),rownames(temp),ignore.case=T),]==0)]
#dim(temp)
#temp<-temp[,which(temp[grep(paste("^","EPCAM","$",sep=""),rownames(temp),ignore.case=T),]==0)]
#dim(temp)

barcode<-colnames(temp)
barcode<-data.frame(Barcode=barcode)
barcode[,1]<-paste(barcode[,1],number,sep="")
write.table(barcode,paste(name,"barcode_filtered.csv",sep="_"),col.names=T,row.names = F,quote = F)
