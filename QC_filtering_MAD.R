################# Di Wu ####################
############## Feb, 2020 ##################

#BiocManager::install("scater")
library("scater")

library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(gplots)
args=commandArgs(TRUE)
path<-args[1]
name<-args[2]
number<-args[3]
species<-args[4]

pbmc.data <- Read10X(data.dir = path)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = name)

##### Using both # of genes  and # of UMIs togather to do filtering ###############
RNA.raw.data <- as.matrix(GetAssayData(pbmc@assays$RNA, slot = "counts"))
  example_sce <- SingleCellExperiment(
    assays = list(counts = RNA.raw.data),
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
  feature.drop <- isOutlier(example_sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
  mt.drop <- isOutlier(test@meta.data$percent.mt, nmads=3, type="higher")
  sce <- example_sce[,!(libsize.drop | feature.drop | mt.drop)]
  
  #sce <- example_sce[,(feature.drop)]
  #sce <- example_sce[,(mt.drop)]
  
  n <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
                  ByMt=sum(mt.drop), Remaining=ncol(sce))
  barcode_MAD <- colnames(sce@assays[["counts"]])
  barcode_MAD<-data.frame(Barcode=barcode_MAD)
  write.table(barcode,paste(name,"_singlet_barcode_filtered_MAD.csv",sep="_"),col.names=T,row.names = F,quote = F)
  NameList_MAD <- as.character(barcode_MAD$Barcode)
  col.num_MAD <- which(colnames(pbmc) %in% NameList_MAD)
  
  NewDF_MAD <- pbmc[,col.num_MAD]
  #NewDF
  #NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^mt-")
  if (species == "human") {
	NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^MT-")
	} else if (species == "mouse") {
		NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^mt-")
	} else if (species == "rat") {
		NewDF_MAD[["percent.mt"]] <- PercentageFeatureSet(NewDF_MAD, pattern = "^Mt-")
	}


  plot1_1_MAD <- FeatureScatter(NewDF_MAD, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2_1_MAD <- FeatureScatter(NewDF_MAD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #saveRDS(NewDF_MAD, paste0("/Users/wud3/Documents/project_tracking/Goodridge_Helen/PB-6925--05--14--2019/merged_run/ADT_HTO/demultiplex/QC_filtering/", name, ".fitered"))
 
  pdf(paste(name,"_QC_scRNA_3MAD.pdf",sep=""),16,12)
  text1_MAD<-paste("Sample Name:",name,sep=" ")
  text2_MAD<-paste(n[1],"cells with a number of UMI counts lower than the median - 3MAD",sep=" ")
  text3_MAD<-paste(n[2],"cells with a number of expressed genes lower than the median - 3MAD.",sep=" ")
  text5_MAD<-paste(n[3],"cells with a % of Mito higher than the median + 3MAD.",sep=" ")
  text4_MAD<-paste("After removing outliers,",n[4],"out of",dim(example_sce)[2],"cells are retained.",sep=" ")
  text_MAD<-paste(text1_MAD,text2_MAD,text3_MAD,text5_MAD,text4_MAD,sep="\n");
  textplot(text_MAD,halign="center",valign="center",cex=2)
  textplot("Before Filtering",halign="center",valign="center",cex=5)
  VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
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
