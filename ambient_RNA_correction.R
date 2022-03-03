#devtools::install_github("constantAmateur/SoupX")
library(SoupX)
library(Seurat)

# Need both raw_feature_bc_matrix and filtered_feature_bc_matrix
args=commandArgs(TRUE)
dataDirs<-args[1]

#Profiling the soup
scl = load10X(dataDirs, keepDroplets = TRUE)
scl$channels$Channel1 = estimateSoup(scl$channels$Channel1)
#modifies the Channel1 SoupChannel object to add estimates of the soup expression profile to each channel entry.
# soup estimation is channel by channel basis

#Visual sanity checks
#data(PBMC_DR)
PBMC_DR <- read.csv("./SoupX_testing/UMAP_50PCs.csv",header = TRUE, row.names = 1)
colnames(PBMC_DR) <- c("RD1", "RD2", "Cluster")
library(ggplot2)
PBMC_DR$GPHA2 = scl$toc["GPHA2", rownames(PBMC_DR)]
gg = ggplot(PBMC_DR, aes(RD1, RD2)) + geom_point(aes(colour = GPHA2 > 0))
plot(gg)
gg = plotMarkerMap(scl, "GPHA2", PBMC_DR)
plot(gg)

#Estimating the contamination fraction
#Picking soup specific genes
scl = inferNonExpressedGenes(scl)

#plot the distribution of expression across cells for the first 20 such genes
tstGenes = rownames(scl$channels$Channel1$nonExpressedGenes)[seq(20)]
gg = plotMarkerDistribution(scl, "Channel1", tstGenes)
plot(gg)

#############################################################################]
# choose appropriate genes for ambient RNA correction #
KRT_genes = c("KRT3", "KRT12", "KRT15", "KRT14","KRT17")
#plot distribution of expression for selected genes
scl = calculateContaminationFraction(scl, "Channel1", list(IG = KRT_genes))
gg = plotChannelContamination(scl, "Channel1")
plot(gg)

library(Matrix)
#rowSums(scl$channels$Channel1$toc[hgGenes, ])
#scl = interpolateCellContamination(scl, "Channel1", useGlobal = TRUE)
scl = interpolateCellContamination(scl, "Channel1", useGlobal = FALSE)
# Need to consider if global correction is appropriate
head(scl$channels$Channel1$rhos)

#Correcting expression profile
#1. produce an expression matrix
scl = strainCells(scl) #creates scl$strainedExp
#2. produce a modified table of counts
#SoupX attempts to remove all the counts that are likely to be soup in origin
scl = adjustCounts(scl) #creates scl$atoc

#Investigating changes in expression
