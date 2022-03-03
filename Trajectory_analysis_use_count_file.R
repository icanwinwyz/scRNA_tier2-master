library(monocle3)
#https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# 3 input files are needed
# 1. expression_matrix, a numeric matrix of expression values, where rows are genes, and columns are cells
# 2. cell_metadata, a data frame, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
# 3. gene_metadata, a data frame, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

#trial data
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

# Make the cell_data_set (CDS) object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#OR
#Generate a cell_data_set from 10X output
# Provide the path to the Cell Ranger output.
#cds <- load_cellranger_data("~/Downloads/10x_data")
#cds <- load_mm_data(mat_path = "~/Downloads/matrix.mtx", 
#                    feature_anno_path = "~/Downloads/features.tsv", 
#                    cell_anno_path = "~/Downloads/barcodes.tsv")


#Pre-process the data
cds <- preprocess_cds(cds, num_dim = 25)

#check that you're using enough PCs to capture most of the variation in gene expression across all the cells in the data set
plot_pc_variance_explained(cds)

#Monocle 3 uses UMAP by default, as we feel that it is both faster and better suited for clustering and trajectory analysis in RNA-seq. 
cds <- reduce_dimension(cds)


#plot_cells(cds)

# color cells by any column of the cell_metadata 
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")

# color cell by expression of interested genes
gene_list = c("che-1",
              "hlh-17",
              "nhr-6",
              "dmd-6",
              "ceh-36",
              "ham-1")
plot_cells(cds, genes=gene_list,label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


#Group cells into clusters
cds = cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
plot_cells(cds, color_cells_by="cell.type")

#Learn the trajectory graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#Order the cells in pseudotime
plot_cells(cds,
          color_cells_by = "embryo.time.bin",
          label_cell_groups=FALSE,
          label_leaves=TRUE,
          label_branch_points=TRUE,
          graph_label_size=1.5)

cds <- order_cells(cds) #will pop up a window and mouse can be used to select root cell
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           label_roots = FALSE)
