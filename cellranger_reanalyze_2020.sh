FOLDER_PATH=$1

cd $FOLDER_PATH/

mkdir $FOLDER_PATH/Final_Results_QC_reanalysis

/usr/bin/Rscript /home/genomics/genomics/bin/10X_scRNA_QC_filtering_Seurat_v3_reanalyze_titan_2020.R $FOLDER_PATH/$2_results/outs/filtered_feature_bc_matrix/ $2 single

mv $FOLDER_PATH/$2_QC_gene300_mt15.pdf $FOLDER_PATH/Final_Results_QC_reanalysis
mv $FOLDER_PATH/$2_QC_MAD.pdf $FOLDER_PATH/Final_Results_QC_reanalysis

/home/genomics/genomics/apps/cellranger-3.1.0/cellranger reanalyze --id=$2_QC_gene300_mt15 --matrix=$2_results/outs/raw_feature_bc_matrix.h5 --barcodes=$2_barcode_filtered_gene300_mt15.csv --nopreflight

/home/genomics/genomics/apps/cellranger-3.1.0/cellranger reanalyze --id=$2_QC_MAD --matrix=$2_results/outs/raw_feature_bc_matrix.h5 --barcodes=$2_barcode_filtered_MAD.csv --nopreflight

cp $FOLDER_PATH/$2_QC_gene300_mt15/outs/*.cloupe $FOLDER_PATH/Final_Results_QC_reanalysis/$2_QC_gene300_mt15.cloupe
cp $FOLDER_PATH/$2_QC_MAD/outs/*.cloupe $FOLDER_PATH/Final_Results_QC_reanalysis/$2_QC_MAD.cloupe
#cp $FOLDER_PATH/$2_results/outs/web* /home/genomics/genomics-archive/NextSeq500_RawData/SampleSheets/Results_Temp/$2_QC.html

mv $FOLDER_PATH/$2_barcode_filtered_gene300_mt15.csv $FOLDER_PATH/Final_Results_QC_reanalysis
mv $FOLDER_PATH/$2_barcode_filtered_MAD.csv $FOLDER_PATH/Final_Results_QC_reanalysis

mv $FOLDER_PATH/$2_Expr_raw_QC_gene300_mt15.csv $FOLDER_PATH/Final_Results_QC_reanalysis
mv $FOLDER_PATH/$2_Expr_raw_QC_MAD.csv $FOLDER_PATH/Final_Results_QC_reanalysis

mv $FOLDER_PATH/$2_Expr_norm_QC_gene300_mt15.csv $FOLDER_PATH/Final_Results_QC_reanalysis
mv $FOLDER_PATH/$2_Expr_norm_QC_MAD.csv $FOLDER_PATH/Final_Results_QC_reanalysis

/home/genomics/bin/change_per.sh $FOLDER_PATH/Final_Results_QC_reanalysis

/home/genomics/bin/change_per.sh $FOLDER_PATH/$2_QC_gene300_mt15
/home/genomics/bin/change_per.sh $FOLDER_PATH/$2_QC_MAD

chmod -R 775 $FOLDER_PATH/Final_Results_QC_reanalysis

chown genomics -R $FOLDER_PATH/Final_Results_QC_reanalysis

echo "Subject: 10X data reanalyze for $2 is done" | sendmail -v di.wu@cshs.org
