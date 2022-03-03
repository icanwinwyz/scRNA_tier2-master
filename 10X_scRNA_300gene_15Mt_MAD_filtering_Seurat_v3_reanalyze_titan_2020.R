###this script is based on Seruat v3.0 and titan server install this package
###### Di Wu #####
###### Mar. 9, 2020 #####
#example: Rscript /home/genomics/genomics/bin/10X_scRNA_QC_filtering_Seurat_v3_reanalyze_titan_2020.R ./Lipo-Fibroblast_results/outs/filtered_feature_bc_matrix/ Lipo-Fibroblast -2
#setwd("~/Documents/Di_scripts/R")
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(gplots)
library(scater)
packageVersion("Seurat")


args=commandArgs(TRUE)
path<-args[1]
name<-args[2]
type<-args[3]

#path<- "./filtered_feature_bc_matrix/"
#name<- "2_test"

test.data <- Read10X(data.dir=path)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 0 cells . Keep all cells with at
# least 300 detected genes
test <- CreateSeuratObject(counts = test.data, min.cells = 0, project = name)
test

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats


test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = paste(c("^MT-", "^Mt-", "^mt-"), collapse="|"))

#human
#test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^MT-")
#Mouse
#test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^mt-")
#rat
#test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^Mt-")

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
pdf(paste(name,"_QC_gene300_mt15.pdf",sep=""),16,12)
if(max(test$percent.mt <= 15 | min(test$nFeature_RNA) >= 300)){
  e <- 0
} else{
  test_4 <- subset(test, subset = nFeature_RNA < 300 & percent.mt > 15)
  e<-dim(test_4@assays$RNA)[2]
}
#test_4 <- subset(test, subset = nFeature_RNA < 300 & percent.mt > 15)
#e<-dim(test_4@assays$RNA)[2]
test_2 <- subset(test, subset = percent.mt <= 15)
a<-dim(test@assays$RNA)[2]-dim(test_2@assays$RNA)[2]-e
test_3 <- subset(test, subset = nFeature_RNA >= 300)
b<-dim(test@assays$RNA)[2]-dim(test_3@assays$RNA)[2]-e
c<-dim(test_1@assays$RNA)[2]
d<-dim(test@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells failed mito% <= 15%",sep=" ")
text3<-paste(b,"cells failed total # expressed genes >= 300.",sep=" ")
text5<-paste(e,"cells failed total # expressed genes >= 300 and mito% <= 15%.",sep=" ")
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


before_QC<-test@assays$RNA

dim(before_QC)

temp<-test_1@assays$RNA
#dim(temp)
#temp<-temp[,which(temp[grep(paste("^","PECAM1","$",sep=""),rownames(temp),ignore.case=T),]==0)]
#dim(temp)
#temp<-temp[,which(temp[grep(paste("^","PTPRC","$",sep=""),rownames(temp),ignore.case=T),]==0)]
#dim(temp)
#temp<-temp[,which(temp[grep(paste("^","EPCAM","$",sep=""),rownames(temp),ignore.case=T),]==0)]
dim(temp)

tag<-NormalizeData(object = test_1)
tag_raw<-GetAssayData(object = tag, slot = "counts")
tag_norm<-GetAssayData(object = tag, slot = "data")
raw_name<-paste(name,"Expr_raw_QC_gene300_mt15.csv",sep="_")
norm_name<-paste(name,"Expr_norm_QC_gene300_mt15.csv",sep="_")

write.csv(tag_raw,raw_name,quote=F,row.names = TRUE)
write.csv(tag_norm,norm_name,quote=F,row.names = TRUE)



#barcode<-colnames(temp)
#barcode<-data.frame(Barcode=barcode)
#barcode[,1]<-paste(barcode[,1],number,sep="")
#write.table(barcode,paste(name,"barcode_filter_single.csv",sep="_"),col.names=F,row.names = F,quote = F)

type <- "single"
if(type == "aggre"){
barcode<-colnames(pbmc_filter@data)
barcode<-data.frame(Barcode=barcode)
write.csv(barcode,paste(name,"barcode_filtered_gene300_mt15.csv",sep="_"),row.names = F,quote = F)
}else if(type == "single"){
barcode<-colnames(temp)
barcode<-data.frame(Barcode=barcode)
barcode[,1]<-paste(barcode[,1],"-1",sep="")
write.csv(barcode,paste(name,"barcode_filtered_gene300_mt15.csv",sep="_"),row.names = F,quote = F)
}

############################## MAD method filtering ##############################

  RNA.raw.data <- as.matrix(GetAssayData(test@assays$RNA, slot = "counts"))
  example_sce <- SingleCellExperiment(
    assays = list(counts = RNA.raw.data)
    #colData = sc_example_cell_info
  )
  
  example_sce <- calculateQCMetrics(example_sce)
  #is.mito <- grepl("^mt-", rownames(example_sce))
  #example_sce <- calculateQCMetrics(example_sce, feature_controls = list(Mt=is.mito))
  #colnames(colData(example_sce))
  #colnames(rowData(example_sce))
  #example_sce
  #plotHighestExprs(example_sce, exprs_values = "counts")
  #plotRowData(example_sce, x = "n_cells_by_counts", y = "mean_counts")
  #plotRowData(example_sce, y="n_cells_by_counts", x="log10_total_counts")
  #plotRowData(example_sce, x="n_cells_by_counts", y="total_counts")
  libsize.drop <- isOutlier(example_sce$total_counts, nmads=3, type="lower", log=TRUE)
  feature.drop <- isOutlier(example_sce$total_features, nmads=3, type="lower", log=TRUE)
  mt.drop <- isOutlier(test@meta.data$percent.mt, nmads=3, type="higher")
  sce <- example_sce[,!(libsize.drop | feature.drop | mt.drop)]
  
  #sce <- example_sce[,(feature.drop)]
  #sce <- example_sce[,(mt.drop)]
  
  n <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
                  ByMt=sum(mt.drop), Remaining=ncol(sce))
  barcode_MAD <- colnames(sce)
  barcode_MAD<-data.frame(Barcode=barcode_MAD)
  barcode_MAD[,1]<-paste(barcode_MAD[,1],"-1",sep="")
  write.table(barcode_MAD,paste(name,"barcode_filtered_MAD.csv",sep="_"),col.names=T,row.names = F,quote = F)
  #NameList_MAD <- as.character(barcode_MAD$Barcode)
  col.num_MAD <- which(paste(colnames(test),"-1",sep="") %in% barcode_MAD$Barcode)

  

  NewDF_MAD <- test[,col.num_MAD]
  #NewDF
  #NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^mt-")
  #if (species == "human") {
#	NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^MT-")
#	} else if (species == "mouse") {
#		NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^mt-")
#	} else if (species == "rat") {
#		NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^Mt-")
#	}
  NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = paste(c("^MT-", "^Mt-", "^mt-"), collapse="|"))

  plot1_1_MAD <- FeatureScatter(NewDF_MAD, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2_1_MAD <- FeatureScatter(NewDF_MAD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #saveRDS(NewDF_MAD, paste0("/Users/wud3/Documents/project_tracking/Goodridge_Helen/PB-6925--05--14--2019/merged_run/ADT_HTO/demultiplex/QC_filtering/", name, ".fitered"))
 
  pdf(paste(name,"_QC_MAD.pdf",sep=""),16,12)
  text1_MAD<-paste("Sample Name:",name,sep=" ")
  text2_MAD<-paste(n[1],"cells with a number of UMI counts lower than the median - 3MAD",sep=" ")
  text3_MAD<-paste(n[2],"cells with a number of expressed genes lower than the median - 3MAD.",sep=" ")
  text5_MAD<-paste(n[3],"cells with a % of Mito higher than the median + 3MAD.",sep=" ")
  text4_MAD<-paste("After removing outliers,",n[4],"out of",dim(example_sce)[2],"cells are retained.",sep=" ")
  text_MAD<-paste(text1_MAD,text2_MAD,text3_MAD,text5_MAD,text4_MAD,sep="\n");
  textplot(text_MAD,halign="center",valign="center",cex=2)
  textplot("Before Filtering",halign="center",valign="center",cex=5)
  VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  CombinePlots(plots = list(plot1, plot2))
  textplot("After Filtering",halign="center",valign="center",cex=5)
  VlnPlot(NewDF_MAD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  CombinePlots(plots = list(plot1_1_MAD, plot2_1_MAD))
  dev.off()
  
NewDF_MAD<-NormalizeData(object = NewDF_MAD)
NewDF_MAD_raw<-GetAssayData(object = NewDF_MAD, slot = "counts")
NewDF_MAD_norm<-GetAssayData(object = NewDF_MAD, slot = "data")
NewDF_MAD_raw_name<-paste(name,"Expr_raw_QC_MAD.csv",sep="_")
NewDF_MAD_norm_name<-paste(name,"Expr_norm_QC_MAD.csv",sep="_")

write.csv(NewDF_MAD_raw,NewDF_MAD_raw_name,quote=F,row.names = TRUE)
write.csv(NewDF_MAD_norm,NewDF_MAD_norm_name,quote=F,row.names = TRUE)
