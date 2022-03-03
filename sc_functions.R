# Useful functions for processing sc data


# QC filtering ------------------------------------------------------------

#' Filter and Summarize Seurat scRNA object
#' 
#' @description 
#' The function takes a Seurat `object` and either determines or applies filtering metrics for mitochondrial percentage, total UMIs, and total genes detected to subset the object. It also saves a .pdf of the cells passing and failing this quality control.
#' 
#' @details 
#' `object` must be a valid Seurat object and should already have a `percent.mt` metadata column summarizing the percent of mitochondrial counts in each cell. This can be calculated with `Seurat::PercentageFeatureSet()`. I am working to extend this to filter a 10X multiome-based seurat object based on relevant parameters, but currently it will only operate on the RNA assay.
#' 
#' If `MAD` is set, all three metrics will have their bounds determined by the mean-absolute deviation method. Note that for some tissue types, this may result in nonsensical cutoffs based on the underlying biology.
#'
#' `thresh.mt`, `thresh.nFeature_RNA`, and `thresh.nCount_RNA` all define the lower and upper bounds of cells that will pass the QC filtering. A typical range for mitochondrial percentages is `c(0,5)`, but this is often violated in many tissues (hearts) or species (humans). Filtering based on detected genes and UMI can be used to remove empty droplets or doublets. Typically 200-300 genes are a minimum. Rangers for UMIs depend on the number of detected genes, but could start around 1000.
#' 
#' The function will output a PDF with several plots summarizing the number of cells passing and failing QC. This file will be written to the current working directory by default, but this can be changed with the `filename` argument. This option along with the `return_type` argument can be useful it you with to apply several different filtering metrics and compare them without actually subsetting the Seurat object that was input.
#' 
#' There is some error handling to control improperly formatter Seurat objects, but this is meant to be an internal function to Cedars-Sinai, so don't go crazy.
#' 
#' @param object a Seurat scRNA-seq object with a `percent.mt` metadata column.
#' @param MAD a boolean value to determine if the mean-absolute deviation should be used to determine cutoffs for filtering.
#' @param thresh.mt a vector of the lower and upper bounds (inclusive) of mitochondrial count percent per cell for cells to retain.
#' @param thresh.nFeature_RNA a vector of the lower and upper bounds (inclusive) of detected genes per cell for cells to retain.
#' @param thresh.nCount_RNA a vector of the lower and upper bounds (inclusive) of UMI count per cell for cells to retain.
#' @param filename a string defining the output filename (with a path if desired)
#' @param return_type should the function return a filtered Seurat object, a list of passing cell barcodes, or just print the summary pdf?
#'
#' @return
#' @export
#'
#' @examples
scRNA_QC <- function(object,
                     MAD = FALSE,
                     thresh.mt = c(0,100),
                     thresh.nFeature_RNA = c(0,Inf),
                     thresh.nCount_RNA = c(0,Inf),
                     filename = "./QC_scRNA_Filtering.pdf",
                     return_type = c("seurat", "csv", "none"))
{
  require(tidyverse)
  require(rcartocolor)
  require(Seurat)
  require(patchwork)
  require(cowplot)
  require(ggvenn)
  # Do some quick error checking
  if(is.null(colnames(object@meta.data))){
    warning("Please specify a valid Seurat object with metadata columns.")
    stop()
  } else if(!"percent.mt" %in% colnames(object@meta.data)) {
    warning("Please run PercentageFeatureSet() for mitochondrial genes and assign result as 'percent.mt'.")
    stop()
  } else if(!"nCount_RNA" %in% colnames(object@meta.data)) {
    warning("nCount_RNA metadata column not found. Are you sure this is a scRNA Seurat object?")
    stop()
  } else if(!"nFeature_RNA" %in% colnames(object@meta.data)) {
    warning("nFeature_RNA metadata column not found. ARe you sure this is a scRNA Seurat object?")
    stop()
  } else if(MAD) {
    message("Filtering thresholds set by mean-absolute deviation. All other filtering thresholds will be ignored.")
  }
  if(MAD==TRUE){
    # Find MAD-based outliers to remove
    require(scater)
    message("Normalizing RNA Assay for MAD-based thresholding.")
    DefaultAssay(object) <- "RNA"
    object <- NormalizeData(object)
    message("Detecting outliers based on MAD.")
    # Remove NA values that sometimes crop up
    object <- object[,!is.na(object$percent.mt)]
    object$filt.mt <- isOutlier(object$percent.mt, 
                                nmads = 3, 
                                type = "higher")
    object$filt.count <- isOutlier(object$nCount_RNA,
                                   nmads = 3,
                                   type = "lower",
                                   log = T)
    object$filt.feat <- isOutlier(object$nFeature_RNA,
                                  nmads = 3,
                                  type = "lower",
                                  log = T)
    # Get MAD thresholds used
    thresh.mt <- c(0, object@meta.data %>% dplyr::filter(filt.mt==FALSE) %>% summarise(max=max(percent.mt), .groups="filt.mt") %>% .$max %>% round(2))
    thresh.nFeature_RNA <- c(object@meta.data %>% dplyr::filter(filt.feat==FALSE) %>% summarise(min=min(nFeature_RNA), .groups="filt.count") %>% .$min, Inf)
    thresh.nCount_RNA <- c(object@meta.data %>% dplyr::filter(filt.count==FALSE) %>% summarise(min=min(nCount_RNA), .groups="filt.feat") %>% .$min, Inf)
  } else if (MAD==FALSE) {
    # Find manually defined outliers to remove
    object$filt.mt <- object$percent.mt < thresh.mt[1] | object$percent.mt >= thresh.mt[2] | is.na(object$percent.mt)
    object$filt.feat <- object$nFeature_RNA < thresh.nFeature_RNA[1] | object$nFeature_RNA >= thresh.nFeature_RNA[2] | is.na(object$nFeature_RNA)
    object$filt.count <- object$nCount_RNA < thresh.nCount_RNA[1] | object$nCount_RNA >= thresh.nCount_RNA[2] | is.na(object$nCount_RNA)
  }
  return_type <- match.arg(return_type)
  # group all cell that are filtered out
  object$filtered <-  object$filt.mt | object$filt.feat | object$filt.count
  # Get cells failing each QC
  MT <- if(sum(object$filt.mt)==0) NULL else Cells(subset(object, subset = filt.mt==TRUE))
  UMIs <- if(sum(object$filt.count)==0) NULL else Cells(subset(object, subset = filt.count==TRUE))
  Genes <- if(sum(object$filt.feat)==0) NULL else Cells(subset(object, subset = filt.feat==TRUE))
  # Plot QC Comparison
  message("Saving filtering plot to ", filename)
  plot_qc <- (
    ( #LHS
      (
        table(object$filtered) %>%
          as.data.frame() %>%
          ggplot(aes(Var1, Freq, col=Var1, fill=Var1)) +
          geom_bar(stat="identity") +
          geom_text(aes(label = Freq), vjust = -0.5, position = position_dodge(0.9)) +
          scale_x_discrete(labels=c("Pass","Fail")) +
          scale_y_continuous(labels = scales::comma,
                             expand = expansion(mult = c(0,.05))) +
          scale_color_carto_d(palette="Safe") +
          scale_fill_carto_d(palette="Safe") +
          theme_cowplot() +
          theme(legend.position = "none",
                plot.title = element_text(hjust=0.5)) +
          labs(title="Passing Cells",y="Cells", x=NULL)
      ) / (
        list(MT=MT,
             UMIs = UMIs,
             Genes = Genes) %>%
          ggvenn(show_percentage = F,
                 fill_color = carto_pal(name="Safe")[c(2,2,2)],
                 stroke_size = 0.5,
                 set_name_size = 4) +
          ggtitle("Failed on") +
          theme(plot.title = element_text(hjust=0.5,
                                          face="bold",
                                          size=12))
      )
    ) | ( #RHS
      (( # Top
        VlnPlot(object,
                features = "nFeature_RNA",
                pt.size = 0,
                group.by = "filtered") +
          labs(x=NULL,y=NULL, title = "Genes")  |
          VlnPlot(object,
                  features = "nCount_RNA",
                  pt.size = 0,
                  group.by = "filtered") +
          labs(x=NULL,y=NULL, title = "UMIs") |
          VlnPlot(object,
                  features = "percent.mt",
                  pt.size = 0,
                  group.by = "filtered") +
          labs(x=NULL,y=NULL, title = "% Mitochondrial")
      ) & scale_fill_carto_d(palette="Safe") &
        scale_y_continuous(labels = scales::comma) &
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank())) / 
        ( # Bottom 
        (
          FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "filtered") + 
            labs(x="UMIs", y="% Mitochondrial") |
            FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "filtered") + 
            labs(x="UMIs", y="Genes") |
            FeatureScatter(object, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "filtered") +
            labs(x="Genes", y="% Mitochondrial")
        ) & theme(plot.title = element_blank()) &
          scale_x_continuous(labels = scales::comma)
      ) & theme(legend.position = "none") & 
        scale_color_carto_d(palette="Safe") &
        scale_y_continuous(labels = scales::comma) 
    )
  ) + plot_layout(widths=c(0.2,1)) +
    plot_annotation(caption = paste0(ifelse(MAD,"MAD-based Filtering Criteria:\n","Manual Filtering Criteria:\n"),
                                     "Genes [",
                                     thresh.nFeature_RNA[1],
                                     ", ",
                                     thresh.nFeature_RNA[2],
                                     "]; UMIs [",
                                     thresh.nCount_RNA[1],
                                     ", ",
                                     thresh.nCount_RNA[2],
                                     "]; % Mitochondrial [",
                                     thresh.mt[1],
                                     ", ",
                                     thresh.mt[2],
                                     "]"))
  print(plot_qc)
  ggsave2(filename,
          height=10,
          width=15)
  if(return_type=="seurat"){
    # output filtered object
    message("Subsetting and returning Seurat object")
    object <- subset(object, filtered==FALSE)
    object$filt.mt <- NULL
    object$filt.feat <- NULL
    object$filt.count <- NULL
    object$filtered <- NULL
    return(object)
  } else if (return_type=="csv"){
    message("Subsetting Seurat object and returning passing cell barcodes")
    object <- subset(object, filtered==FALSE)
    return(Cells(object))
  }
}



# cluster plots -------------------------------------------------------------------
# to generate the typical tsne and umap plots with and without cluster labels
viz_clusters <- function(object,
                         group.var = "seurat_clusters",
                         group.label = "Cluster",
                         path = "./",
                         label.plot=c("both","labeled", "unlabeled"),
                         Visium=FALSE,
                         cols = NULL,
                         ...)
{
  # load require packages
  label.plot <- match.arg(label.plot)
  require(ggplot2)
  require(tidyverse)
  require(rcartocolor)
  require(Seurat)
  require(patchwork)
  require(cowplot)  
  # Check if it's a visium object and plot the slide images too
  invisible(if(Visium){
    message("Visium object found. Plotting ", group.var, " on slide images.")
    pv <- SpatialDimPlot(object,
                         group.by = group.var) %>% wrap_plots(guides = "collect") &
      {if(!is.null(cols))scale_color_manual(values=cols) else if(dim(unique(object[[group.var]]))[1]<13) scale_color_carto_d(name=group.label,palette = "Safe")}
    print(pv)
    ggsave2(filename = paste0(path,
                              "Slides_by_",
                              group.label,
                              ".pdf"),
            ...)
  })
  invisible(lapply(c("tSNE", "UMAP"),
                   function(p){
                     # compensate for long labels
                     width.adj <- 6+0.077*max(nchar(as.character(unlist(unique(object[[group.var]])))))
                     if(label.plot!="unlabeled"){
                       #Plot labeled
                       p1 <- DimPlot(object,
                                     group.by = group.var,
                                     label=T,
                                     repel=T,
                                     reduction = tolower(p)) +
                         labs(x=paste(p,"1"),
                              y=paste(p,"2"),
                              title=paste(p,"by", group.label)) +
                         {if(!is.null(cols))scale_color_manual(values=cols) else if(dim(unique(object[[group.var]]))[1]<13) scale_color_carto_d(name=group.label,palette = "Safe")}  +
                         scale_x_continuous(breaks = NULL) +
                         scale_y_continuous(breaks = NULL)
                       print(p1)
                       ggsave2(filename = paste0(path,
                                                 p,
                                                 "_by_",
                                                 group.label,
                                                 ifelse(label.plot!="both",
                                                        "",
                                                        "_Labeled"),
                                                 ".pdf"),
                               height=6,
                               width=width.adj)
                     }
                     #Plot unlabeled
                     if(label.plot!="labeled"){
                       p2 <- DimPlot(object,
                                     group.by = group.var,
                                     label=F,
                                     reduction = tolower(p)) +
                         labs(x=paste(p,"1"),
                              y=paste(p,"2"),
                              title=paste(p,"by", group.label)) +
                         {if(!is.null(cols))scale_color_manual(values=cols) else if(dim(unique(object[[group.var]]))[1]<13) scale_color_carto_d(name=group.label,palette = "Safe")}  +
                         scale_x_continuous(breaks = NULL) +
                         scale_y_continuous(breaks = NULL)
                       print(p2)
                       ggsave2(filename = paste0(path,
                                                 p,
                                                 "_by_",
                                                 group.label,
                                                 ifelse(label.plot!="both",
                                                        "",
                                                        "_Unlabeled"),
                                                 ".pdf"),
                               height=6,
                               width=width.adj)
                     }
                   }))
}


