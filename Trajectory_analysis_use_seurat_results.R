library(xlsx)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(monocle3)

########## Extract cells from datasets ######################################################
load("E9_5_11_5.rds")
DimPlot(E9_5_11_5, reduction = "umap", group.by = "group", pt.size = 1.5)
DimPlot(E9_5_11_5, pt.size = 2.0) # clusters
dim(E9_5_11_5[["RNA"]]@data)
# 13285   215
E9_5_11_5 = subset(E9_5_11_5,idents = 'endoderm')
table(E9_5_11_5@meta.data$group)
# Blanca_lung_E10_5    Blanca_lung_E11_5 Blanca_lung_E9_5_em2 Blanca_lung_E9_5_em3 
#      22                   22                   8                   14
E9_5_11_5 = RenameIdents(E9_5_11_5, 'endoderm' = 'E9.5_11.5_Endoderm')


load("Cohen_E12_5_CD45N_Mes.rds")
DimPlot(Cohen_E12_5_CD45N_Mes, reduction = "umap", group.by = "group", pt.size = 1.5)
DimPlot(Cohen_E12_5_CD45N_Mes, pt.size = 2.0) # clusters
dim(Cohen_E12_5_CD45N_Mes[["RNA"]]@data)
# 14136  1158
Cohen_E12_5_CD45N_Mes = subset(Cohen_E12_5_CD45N_Mes,idents = c('Pro-lipo','Myo','Inter'))
table(Cohen_E12_5_CD45N_Mes@meta.data$group)
# Cohen_E12_5_CD45N1 Cohen_E12_5_CD45N2 Cohen_E12_5_CD45N3 Cohen_E12_5_CD45N4 
#       266                258                263                279
Cohen_E12_5_CD45N_Mes = RenameIdents(Cohen_E12_5_CD45N_Mes, 'Pro-lipo' = 'E12.5_Pro-lipo','Myo' = 'E12.5_Myo','Inter'='E12.5_Inter')

  
load("E17_5_Wt_Fib_cells.rds")
DimPlot(Fib.cells, reduction = "umap", group.by = "group", pt.size = 1.5)
DimPlot(Fib.cells, pt.size = 2.0) # clusters
dim(Fib.cells[["RNA"]]@data)
# 16315  2002
Fib.cells = subset(Fib.cells,idents = c('Lipo','Myo','Ebf1+','Inter'))
table(Fib.cells@meta.data$group)
#  WT1 WT2 WT3 
#  705 565 581
Fib.cells = RenameIdents(Fib.cells, 'Lipo' = 'E17.5_Lipo','Myo' = 'E17.5_Myo','Ebf1+'='E17.5_Ebf1+','Inter'='E17.5_Inter')

######### Integrate Data using Seurat #################################
All_data.anchors <- FindIntegrationAnchors(object.list = list(E9_5_11_5,Cohen_E12_5_CD45N_Mes,Fib.cells), dims = 1:20) #anchor.features = 2000
All_data <- IntegrateData(anchorset = All_data.anchors, dims = 1:20)
#All_data <- IntegrateData(anchorset = All_data.anchors, dims = 1:30, normalization.method = "SCT",new.assay.name = "integrated_SCT")

DefaultAssay(All_data) <- "integrated"
All_data <- ScaleData(All_data, verbose = FALSE)
All_data <- RunPCA(All_data, npcs = 30, verbose = FALSE)
#All_data <- RunTSNE(All_data, npcs = 30, verbose = FALSE)
All_data <- RunUMAP(All_data, reduction = "pca", dims = 1:20, min.dist = 0.3) #min.dist = 1.0
#All_data <- RunUMAP(All_data, reduction = "pca", dims = 1:20, min.dist = 1.0, n.components = 3)

png("UMAP_Groups_Seurat_Integration.png",width = 7, height = 5, units = 'in', res = 300)
DimPlot(All_data, reduction = "umap", group.by = "group", pt.size = 1.0)
dev.off()
#DimPlot(All_data, reduction = "tsne", group.by = "group", pt.size = 1.0)
#DimPlot(All_data, reduction = "pca", pt.size = 1.0) # group.by = "group"



Idents(All_data) = factor(Idents(All_data),levels = c("E9.5_11.5_Endoderm","E12.5_Pro-lipo","E17.5_Lipo",
                                                      "E12.5_Myo","E17.5_Myo","E12.5_Inter","E17.5_Inter",
                                                      "E17.5_Ebf1+"))
levels(Idents(All_data))





png("UMAP_Clusters_Seurat_Integration.png",width = 7, height = 5, units = 'in', res = 300)
DimPlot(All_data,reduction = "umap", pt.size = 1.0) # clusters
#DimPlot(All_data,reduction = "tsne", pt.size = 1.0) # clusters
dev.off()

save(All_data, file = "E9.5_E11.5_E12.5_E17.5.rds")
load("E9.5_E11.5_E12.5_E17.5.rds")

################## run monocle3 for trajectories ####################################################
#### convert Seurat 3 object to Monocle 3 object
### Building the necessary parts for a basic cds

seurat = All_data
# part one, gene annotations

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix

New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
#cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- Idents(seurat)
  
  
### Assign the cluster info

list_cluster <- Idents(seurat)
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

### Could be a space-holder, but essentially fills out louvain parameters

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

### Assign UMAP coordinate

cds_from_seurat@reducedDims@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings


### Assign feature loading for downstream module analysis

cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

############# add clusters info. 
colData(cds_from_seurat)$clusters = list_cluster

### Learn graph, this step usually takes a significant period of time for larger samples

print("Learning graph, which can take a while depends on the sample")

cds_from_seurat <- learn_graph(cds_from_seurat,use_partition = T) # use_partition = T
#plot_cells(cds_from_seurat, label_groups_by_cluster=F,  color_cells_by = "cluster",group_label_size = 4,graph_label_size = 4,
#           cell_size = 1, trajectory_graph_segment_size = 1,label_leaves=FALSE,label_branch_points=FALSE)

##### order cells #######################################
#root_cell_list <- grep("lung_E9.5_embryo3", counts(cds_from_seurat)@Dimnames[[2]])
#root_cell_list <- counts(cds_from_seurat)@Dimnames[[2]][root_cell_list]

cds_from_seurat <- order_cells(cds_from_seurat, root_cells = "lung_E9.5_embryo3_sc3") # root_cell_list

png("UMAP_Clusters_Pseudotime_Seurat_Integration.png",width = 7, height = 5, units = 'in', res = 300)
plot_cells(cds_from_seurat, label_groups_by_cluster=F,  color_cells_by = "pseudotime",label_cell_groups=FALSE,
           group_label_size = 4,graph_label_size = 4,
           cell_size = 1, trajectory_graph_segment_size = 1.2,label_leaves=FALSE,label_branch_points=F)
dev.off()
png("UMAP_Clusters_Monocle3_Seurat_Integration.png",width = 8, height = 5, units = 'in', res = 300)
plot_cells(cds_from_seurat, label_groups_by_cluster=F,  color_cells_by = "cluster",label_cell_groups=FALSE,
           group_label_size = 4,graph_label_size = 4,
           cell_size = 1, trajectory_graph_segment_size = 1.2,label_leaves=FALSE,label_branch_points=F)
dev.off()
