# scRNA_tier2  

Trial data can be donwloaded at    
https://genomicscloud.csmc.edu/apps/files/?dir=%2FGenomicsOwnCloud%2FMarban_Eduardo%2FGDC-6864--05--06--2019%2FCtrl-1_results%2Ffor_secondary_analysis%2Ffiltered_feature_bc_matrix

### Generating two sets of QC results using MAD and gene=<300 and/or mito%>=15% (on Titan)  ##
```
source activate my_python_3_7

Rscript /home/genomics/genomics/bin/10X_scRNA_300gene_15Mt_MAD_filtering_Seurat_v3_reanalyze_titan_2020.R /path/to/filtered_feature_bc_matrix sample1

eg: Rscript /home/genomics/genomics/bin/10X_scRNA_QC_filtering_Seurat_v3_reanalyze_titan_2020.R /home/genomics/genomics/data/Temp/Sequence_Temp/NovaSeq/Fastq_Generation/200429_A00319_0145_AHLYF7DRXX_11_10_22/BD-9430--04--14--2020-_Knott_Simon_scRNA_10X_Human/JasBioBankProstatetumorpunch_results/outs/filtered_feature_bc_matrix/ JasBioBankProstatetumorpunch

```

![Workflow](https://github.com/cedars-sinai-genomics-core/scRNA_tier2/blob/master/scRNA-seq_tier2_workflow.png)

The reanalysis step can also run through Titan portal: 
http://10.220.239.17/upload_cellranger_reanalyze.php

Papers using MAD method for QC filtering:

https://www.sciencedirect.com/science/article/pii/S2211124718320746 (Cell Reports)  
https://www.nature.com/articles/s41591-019-0468-5 (Nat Med)  
https://www.nature.com/articles/s41593-019-0393-4 (Nat Neurosci)  

   
