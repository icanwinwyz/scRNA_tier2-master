FROM satijalab/seurat:4.0.3

# Install CRAN packages alphabetically
RUN install2.r --error \
	ggvenn \
	rcartocolor \
	SoupX \
	UpSetR \
	tidyverse 

# Install Bioconductor packages in order of need
RUN R -e 'BiocManager::install(c("scater", "biomaRt", "SingleR", "celldex", "DOSE", "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db", "org.Cf.eg.db", "DropletUtils"), ask=FALSE, update=TRUE)'

RUN apt-get update && apt-get -y install python3-pip

RUN pip install kb-python cellrank
