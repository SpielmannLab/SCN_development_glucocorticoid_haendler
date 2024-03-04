#!/usr/bin/bash

file_sc_obj="/data/humangen_singlecell/Astiz_SCN_development/figures/clusters/9_seurat_dim_red/merged.rds_dim_reduction.rds"
outdir="/data/humangen_mouse/2023_glucocorticoid/output/figure1_expression_featureplot/"
mkdir -p $outdir

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

conda activate AstizSCN

./figure1_expression_featureplot.R $file_sc_obj "Syt1" "white" "blue" $outdir &
./figure1_expression_featureplot.R $file_sc_obj "Pdgfra" "white" "green" $outdir &
./figure1_expression_featureplot.R $file_sc_obj "Aldh1l1" "white" "red" $outdir &

./figure1_expression_featureplot.R $file_sc_obj "Syt1" "black" "blue" $outdir &
./figure1_expression_featureplot.R $file_sc_obj "Pdgfra" "black" "green" $outdir &
./figure1_expression_featureplot.R $file_sc_obj "Aldh1l1" "black" "red" $outdir 