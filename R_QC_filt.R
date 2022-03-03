# Example usage
# positional arguments
# 1: path to cellranger filtered_feature_bc_natrix.h5 file
# 2: whether to apply MAD based filtering or not
# 3: filename (with path) for the output QC plot
# 4: filename (with path) for the filtered Seurat RDS
# 5-7: if $2 is FALSE, mt, feature, and UMI ranges (in order) for cells to retain
#
# This will take a cellranger output, apply manual cutoffs of 0-6% for mt reads, 300-Inf featurs and 1000-Inf UMIs, save a list of passing barcodes, and output the summary plot
# Rscript ./R_QC_filt.R ~/path/to/10x/filtered_feature_bc_matrix.h5 FALSE ./Plot.pdf ./Barcodes.csv 0-6 300-Inf 1000-Inf 
#
# This will do the same as above except with MAD filtering
# Rscript ./R_QC_filt.R ~/path/to/10x/filtered_feature_bc_matrix/ TRUE ./Plot.pdf ./Barcodes.csv  


source("./sc_functions.R")

suppressPackageStartupMessages({
  require(tidyverse)
  require(rcartocolor)
  require(Seurat)
  require(patchwork)
  require(cowplot)
  require(ggvenn)
})


# Error checking for command arguments
args <- commandArgs(trailingOnly = T)

if (length(args)<4) {
  stop("Please specify the necessary command arguments. Check the header of this script for more info.",
       call.=FALSE)
}
## check input file
if(file.exists(args[1])) {
  path.h5 <- args[1]
} else {
  stop("Input h5 file does not exist.",
       call. = FALSE)
}
## check MADS argument
if(!is.na(as.logical(args[2]))) {
  bool.MAD <- as.logical(args[2])
} else {
  stop("Specify TRUE/FALSE for MAD filtering.",
       call. = FALSE)
}
## are outputs writeable?
if ((file.access(dirname(args[3]))*file.access(dirname(args[4])))==0) {
  path.plot <- args[3]
  path.csv <- args[4]
} else {
  stop("One or both outputs are not in writable directories.",
       call. = FALSE)
}
## manual thresholds if necessary
if (bool.MAD & length(args)>4) {
  stop("You must set either MAD filtering to TRUE or supply all manual thresholds. These options are mutally exclusive.")
} else if (!bool.MAD & length(args)>4 & length(args)<7) {
  stop("Currently you must specify criteria for all three manual thresholds.")
} else if (!bool.MAD & length(str_split(args[5:7], pattern = "-", simplify = T))>6) {
  stop("Manual threshold ranges must separated by hyphens, e.g. 0-7 200-Inf 1000-Inf")
} else if (!bool.MAD & any(is.na(as.numeric(str_split(args[5:7], pattern = "-", simplify = T))))) {
  stop("Manual threshold ranges can only contain numeric values or Inf")
} else {
  thresholds <- as.numeric(str_split(args[5:7], pattern = "-", simplify = T))
}

# read in and create the seurat object
message("Creating Seurat object from ", basename(args[1]))
counts <- Read10X_h5(args[1])
object <- CreateSeuratObject(counts=counts,
                             project = "scRNAseq",
                             min.cells = 3) %>% 
  PercentageFeatureSet(pattern = "^[Mm][Tt]-", col.name = "percent.mt")

message("Filtering and saving summary plot to ", path.plot)
object <- scRNA_QC(object,
                   MAD = bool.MAD,
                   thresh.mt = thresholds[c(1,4)],
                   thresh.nFeature_RNA = thresholds[c(2,5)],
                   thresh.nCount_RNA = thresholds[c(3,6)],
                   filename = path.plot,
                   return_type = "csv")

message("Saving passing cell barcodes to ", path.csv)
write.table(object,
            row.names = F,
            col.names = F,
            file = path.csv)

# remove the Rplots.pdf that is made by mistake
invisible(if(file.exists("./Rplots.pdf")){
  invisible(file.remove("./Rplots.pdf"))
})





